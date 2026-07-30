%% DRM indexing engine - pure alpha-Ti (HCP)  ------------------------------
% Parallel pipeline to main.m, rebuilt for hexagonal (6/mmm) crystallography.
% New / HCP-specific functions:
%   facetNormalsHCP  - {hkil} families -> correct crystal-frame facet normals
%   rotate_facet_hcp - rotate those normals into the specimen frame
%   DRPsim_hcp       - simulate a DRP from HCP facet normals
%   makeDRPdic_hcp   - orientation dictionary (hexagonal fundamental zone)
%   DRPLibGenerator_hcp / check_indexing_result_hcp
% -------------------------------------------------------------------------
addpath(genpath(fileparts(mfilename('fullpath'))));

addpaths=false;
if addpaths
    addpath("C:\Users\mrbla\Desktop\Cambridge\Physics-based DRM\DRM_MATLAB\baseFunctions") % add path to base functions
    addpath("C:\Users\mrbla\Desktop\Cambridge\Physics-based DRM\DRM_MATLAB\mathFunctions")
    addpath("C:\Users\mrbla\Desktop\Cambridge\Physics-based DRM\DRM_MATLAB\utilities")
    addpath("C:\Users\mrbla\Desktop\mrbla-downloads\Relevant Downloads\MTEX Code\mtex-6.1.0") % for MTEX
    startup_mtex
end

%% experiment (illumination) parameters -----------------------------------
exp_para.th_max = 65;
exp_para.th_min = 10;
exp_para.th_num = 12;
exp_para.ph_num = 120;
exp_para.ph_min = 0;
exp_para.ph_max = 357;

% peak / band shape:  [i_Main, i_facet, sd_Main, sd_facet]
%   sd_Main = peak half-width (deg), sd_facet = band half-width (deg)
exp_para.fitting_para = [1, 0.6, 7, 6];
exp_para.peakModel = "specular";   % "specular" (sharp, tunable) or "cosine"

pos1 = [0 0 0 0];                   % ROI; [0 0 0 0] -> draw one interactively
scaleCoeff   = 0.5;
use_saved_drp_dic  = false;
use_autoencoder    = false;

%% crystallography ---------------------------------------------------------
% One crystalSymmetry object drives faceting, the orientation grid AND the IPF
% colouring, so every stage is on the same convention.
exp_para.crystal = crystalSymmetry('6/mmm',[2.95 2.95 4.68], ...
    'X||a*','Y||b','mineral','Ti-alpha');

% --- pick the facet family(-ies) that etch on YOUR sample --------------
% Each row is a Miller-Bravais plane {h k i l}, i = -(h+k).  List several rows
% to model co-existing families; set weights to bias their peak/band strength.
% Determine the right family by comparing a grain of KNOWN orientation (from
% EBSD) against DRPsim_hcp for each candidate below.
exp_para.facetHKIL = [ 1  0 -1  0 ];   % first-order prism {10-10}  (default)
% exp_para.facetHKIL = [ 1  1 -2  0 ];   % second-order prism {11-20}
% exp_para.facetHKIL = [ 1  0 -1  1 ];   % first-order pyramidal {10-11}
% exp_para.facetHKIL = [ 0  0  0  1 ];   % basal (0001)
% exp_para.facetHKIL = [ 1 0 -1 0; 1 0 -1 1 ];   % prism + pyramidal combined
% exp_para.facetFaceWeights = [1; 0.5];          % (optional) per-family peak weights
% exp_para.facetBandWeights = [1; 0.5];          % (optional) per-family band weights

exp_para = facetNormalsHCP(exp_para);            % builds exp_para.facetNormals etc.

% --- surface-topology (reflector) model --------------------------------
% AFM/SEM show smooth elongated RIDGES, not flat facets.  Ridges are the SURFACE
% TRACES of a crystal plane (axis t = m x Z): t is horizontal, so the DRP stripe
% ALWAYS passes through the pattern centre (as observed).  Model:
%   "facet" flat facets (peaks+bands) | "ridge" cylinders (stripes) | "both"
exp_para.reflectorModel = "ridge";     % <- "facet" | "ridge" | "both"
exp_para.ridgeMode      = "planeTrace";% "planeTrace" (trace m x Z) | "direction" old
exp_para.ridgeFitting   = [1, 6];      % [stripe amplitude, cross-stripe width deg]
exp_para.ridgeTaper     = 0;           % >0 fades the stripe toward the DRP edge
% ridge PLANE whose surface traces form the ridges; pick via compareRidgeAxes.
% A single plane (e.g. basal) -> a single stripe through the centre.
exp_para.ridgeHKIL = [ 0  0  0  1 ];   % basal (0001)  (default: one stripe)
% exp_para.ridgeHKIL = [ 1  0 -1  0 ];   % prism {10-10}
% exp_para.ridgeHKIL = [ 1  0 -1  1 ];   % pyramidal {10-11}
% exp_para.ridgeHKIL = [ 1  0 -1  2 ];   % pyramidal {10-12}
exp_para = ridgePlanesHCP(exp_para);   % builds exp_para.ridgePlaneNormals etc.
%
% Alternative "direction" geometry (rigid crystal direction, off-centre stripe):
%   exp_para.ridgeMode = "direction";
%   exp_para.ridgeUVTW = [1 1 -2 0];  exp_para = ridgeAxesHCP(exp_para);

%% load sample dataset -----------------------------------------------------
[igray_sample, phitheta, pos, img_sample] = drp_loader( ...
    exp_para, pos1, format='jpg', scale=scaleCoeff);
igray_norm = igray_sample;

% convert into DRP stack  [n1 x n2] cells, each th_num x ph_num
drp_original = igray2drp(igray_norm, phitheta, exp_para);

%% (do this FIRST) A/B the facet families on a known grain -----------------
% Pick exp_para.facetHKIL by eye before committing to a dictionary.  Give the
% EBSD orientation(s) of one or more grains and (ideally) their measured DRP,
% and compare candidate families side by side.

% d1 till d5
euKnown = [26.7 55.5 20.6; 27.9 67.5 19.0; 175.6 122.0 59.5; 37.0 84.7 31.9; 0.1 42.0 1.5; 11.4 50.7 2.3; 40.9 75.8 1.3];                       % [phi1 Phi phi2] from EBSD (deg)
%euKnown = [26.7 55.5 20.6; 27.9 67.5 19.0; 11.4 50.7 2.3; 40.9 75.8 1.3]; 
%measDRP = {drp_original{359,110}, drp_original{350, 330}, ...
% drp_original{320, 370}, drp_original{350, 370}, drp_original{620, 345}};           % measured DRP of that same grain - row col
d1 = drp_original{110, 359}; % reverse order of DRP viewer
d2 = drp_original{330, 350};
d3 = drp_original{370, 320};
d4 = drp_original{370, 350};
d5 = drp_original{345, 620};
d6 = drp_original{405, 109};
d7 = drp_original{490, 1110};
%d5 = drp_original{};
%compareFacetFamilies(euKnown, exp_para, measured={measDRP});
compareFacetFamilies(euKnown, exp_para, measured={d1;d2;d3;d4; d5; d6; d7});

% Multiple grains at once (rows = grains), default 4 candidate families:
%   euKnown = [30 45 0; 10 20 0; 0 80 30];
%   compareFacetFamilies(euKnown, exp_para, measured={d1; d2; d3});
%
% Custom candidate list (e.g. prism vs pyramidal vs prism+pyramidal):
%   fam = { [1 0 -1 0], [1 0 -1 1], [1 0 -1 0; 1 0 -1 1] };
%   compareFacetFamilies(euKnown, exp_para, families=fam);
%
% Then set exp_para.facetHKIL above to the winner and re-run facetNormalsHCP.

%% (ridge theory) A/B candidate ridge planes on the known grain ------------
% Tests the ridge -> stripe model (planeTrace): which crystal PLANE's surface
% trace reproduces the measured streak?  Columns = candidate {hkil}; every
% stripe passes through the centre (t = m x Z is horizontal), as observed.
%compareRidgeAxes(euKnown, exp_para, measured={measDRP});
compareRidgeAxes(euKnown, exp_para, measured={d1;d2;d3; d4; d5; d6; d7});
%   compareRidgeAxes(euKnown, exp_para, measured={measDRP}, mode="direction");  % old off-centre geom

% Confirm the plane + EBSD->DRM registration: overlay the predicted basal stripe
% on each measured DRP, and sweep phiOffset until ONE value aligns ALL grains.
checkRidgeAzimuth(euKnown, {d1;d2;d3;d4; d5; d6; d7}, exp_para, ridgeHKIL=[0 0 0 1], phiOffset=0);
% e.g. try phiOffset = 0:30:150 to find the global mounting offset:
for off = 0:30:150
    checkRidgeAzimuth(euKnown, {d1;d2;d3;d4; d5; d6; d7}, exp_para, phiOffset=off);
end

%% (multi-plane theory) superimpose several planes' DRPs for a grain -------
% Test whether basal + another plane reproduces each measured DRP.  Each entry
% = {hkil, model, variant, weight}: model "ridge"/"facet", variant "family"
% (all symmetry equivalents) or "single" (that exact plane -> variant selection).
planes = { {[0 0 0 1], "ridge", "family", 1.0}, ...    % basal ridge
           {[1 0 -1 0], "ridge", "family", 0.7}, ...   % prism ridges
           {[1 0 -1 1], "ridge", "family", 0.5} };     % pyramidal ridges
superimposePlaneDRPs(euKnown, exp_para, planes, measured={d1;d2;d3;d4;d5;d6;d7});
% incoherent addition instead of brightest-wins:
%   superimposePlaneDRPs(euKnown, exp_para, planes, measured={d1;...}, combine="sum");
% variant selection - one SPECIFIC prism plane rather than the whole family:
%   planes = { {[0 0 0 1],"ridge","family",1}, {[1 0 -1 0],"ridge","single",0.8} };

18181818

% Sharpen (narrower stripe) and taper it toward the centre to match a streak:
%   compareRidgeAxes(euKnown, exp_para, measured={measDRP}, ...
%                    ridgeFitting=[1 4], ridgeTaper=6);
%
% Facet vs ridge vs both for one orientation (direct visual):
%   ep_f = exp_para; ep_f.reflectorModel = "facet";
%   ep_r = exp_para; ep_r.reflectorModel = "ridge";
%   ep_b = exp_para; ep_b.reflectorModel = "both";
%   figure, tiledlayout(1,4,'TileSpacing','compact')
%   nexttile, DRPdisp(DRP_norm(measDRP),exp_para);                                   title('measured')
%   nexttile, DRPdisp(DRPsim_hcp(euKnown(1),euKnown(2),euKnown(3),ep_f),exp_para);   title('facet')
%   nexttile, DRPdisp(DRPsim_hcp(euKnown(1),euKnown(2),euKnown(3),ep_r),exp_para);   title('ridge')
%   nexttile, DRPdisp(DRPsim_hcp(euKnown(1),euKnown(2),euKnown(3),ep_b),exp_para);   title('both')
%
% Happy with the ridge model? Set exp_para.reflectorModel above to "ridge" (or
% "both") and the dictionary/indexing below use it automatically.

%% DRP dictionary ----------------------------------------------------------
if ~use_saved_drp_dic
    [drpDic, euDic, rotDic] = makeDRPdic_hcp(exp_para, resolution=3);
    save('DRP_dictionary_PureTi_HCP.mat','drpDic','euDic','rotDic');
else
    S = load('DRP_dictionary_PureTi_HCP.mat','drpDic','euDic','rotDic');
    drpDic = S.drpDic; euDic = S.euDic; rotDic = S.rotDic;
end

%% autoencoder -------------------------------------------------------------
if ~use_autoencoder
    AE_DRM = trainAutoencoder(drpDic, 100, ...
        'MaxEpochs',200, 'L2WeightRegularization',0.001, ...
        'SparsityRegularization',4, 'SparsityProportion',0.10, ...
        'ScaleData',false, 'UseGPU',false);
    save('DRP_autoencoder_PureTi_HCP.mat','AE_DRM');
else
    S = load('DRP_autoencoder_PureTi_HCP.mat','AE_DRM');
    AE_DRM = S.AE_DRM;
end

%% index -------------------------------------------------------------------
tic
index_result = IndexingEngine(drp_original, AE_DRM, exp_para, drpDic, euDic, rotDic);
toc

%% results -----------------------------------------------------------------
[n1,n2] = size(drp_original);
index_num = cellfun(@(x) sum(x,'all'), drp_original);
non_index_bg = index_num > 3e4;
figure, imshow(non_index_bg); title('background mask (tune threshold)');

figure, imshow(plot_ipf_map(index_result.EUmap)); title('IPF-z, pure Ti (HCP)');

%% validate: measured vs predicted DRP at clicked pixels -------------------
[drp_measurement, drp_predicted, xy] = check_indexing_result_hcp( ...
    index_result.EUmap, drp_original, exp_para);

%% (optional) direct indexing without autoencoder -------------------------
% drpLib = DRPLibGenerator_hcp(5*degree, exp_para);
% indexResult = DirectDIEngine(drp_original, drpLib);
% figure, imshow(plot_ipf_map(indexResult.euMap));
% title('direct indexing (no autoencoder), pure Ti (HCP)');
