%% Physics-based DRM indexing - WROUGHT Ti-6Al-4V, alpha (HCP) phase --------
%
% Implements the etching / DRP-generation mechanism documented in section
% 2.3.2 of the thesis chapter ("Directional reflectance patterns", Figs 2.11 -
% 2.15) rather than the flat-facet model used for cubic alloys or the single
% basal-trace model used for pure alpha-Ti (main_PureTi_HCP.m).
%
% THE MECHANISM IN ONE PARAGRAPH ------------------------------------------
% Wrought Ti64 is a stack of ~2 um alpha laths separated by a ~150 nm fine
% beta matrix.  Kroll's reagent etches alpha much faster than beta, and within
% alpha the basal plane etches slowest.  The etched surface is therefore a
% lamellar corrugation whose crests are ROUNDED (semi-cylindrical), so each
% crest behaves as a convex mirror: it fans the incoming light out along the
% great circle perpendicular to its axis, which appears in the DRP as a STRIPE
% through the pattern centre.  Two crest families exist, at right angles:
%
%   * the ridge joining the exposed PRISM facet and the broad face of a lath runs PARALLEL
%     to the c-axis (all prism planes contain c), so it makes a stripe
%     PERPENDICULAR to the c-axis direction.  It needs the prism plane and broad face to be
%     exposed -> switches on at HIGH Phi.  Dominant at Fig 2.12 point 1
%     (Phi = 87.9 deg).
%
%   * the rounded edge between the beta matrix and the slow-etching BASAL
%     plane is PERPENDICULAR to the c-axis, so
%     it makes a stripe PARALLEL to the c-axis direction.  It needs basal and broad face
%        to
%     be exposed -> switches on at LOW Phi.  Dominant at Fig 2.12 point 4
%     (Phi = 37.8 deg).
%
% Both are always present; only their intensity ratio changes with Phi and phi2, which
% is exactly the observation that DRPs show two stripes of unequal magnitude.
% Two point-like features complete the picture: the broad plane and basal plane acts as
% a plane mirror,
%
% WHERE THE CODE LIVES ----------------------------------------------------
%   lathGeometryTi64            model setup: crystallography + all parameters
%   lathReflectors              THE PHYSICS, per orientation, fully exposed
%   DRPsim_Ti64                 paints the reflectors into a DRP
%   checkTi64Mechanism          verification: decompose + overlay + print
%   makeDRPdic_Ti64             orientation dictionary
%   check_indexing_result_Ti64  measured vs predicted at clicked pixels
%   borSiblingRotations         Burgers-OR sibling variants (optional)
%
% Read lathReflectors.m top to bottom to audit the mechanism; every
% intermediate quantity (tilt, exposure, coincidence, stripe azimuth,
% specular theta) is returned and printed rather than buried.
% -------------------------------------------------------------------------
addpath(genpath(fileparts(mfilename('fullpath'))));

addpaths = false;
if addpaths
    addpath("C:\Users\mrbla\Desktop\Cambridge\Physics-based DRM\DRM_MATLAB\baseFunctions")
    addpath("C:\Users\mrbla\Desktop\Cambridge\Physics-based DRM\DRM_MATLAB\mathFunctions")
    addpath("C:\Users\mrbla\Desktop\Cambridge\Physics-based DRM\DRM_MATLAB\utilities")
    addpath("C:\Users\mrbla\Desktop\mrbla-downloads\Relevant Downloads\MTEX Code\mtex-6.1.0")
    startup_mtex
end

%% experiment (illumination) parameters ------------------------------------
% From section 2.3.2: 3 deg in phi, 5 deg in theta, theta in [10,65]
% -> 120 x 12 = 1440 images.
exp_para.th_min = 10;
exp_para.th_max = 65;
exp_para.th_num = 12;
exp_para.ph_min = 0;
exp_para.ph_max = 357;
exp_para.ph_num = 120;

% kept for compatibility with the shared helpers (unused by DRPsim_Ti64,
% which takes its widths from exp_para.lath instead)
exp_para.fitting_para = [1, 0.6, 7, 6];

pos1              = [0 0 0 0];    % ROI; [0 0 0 0] -> draw one interactively
scaleCoeff        = 0.5;
use_saved_drp_dic = false;
use_saved_ae      = false;

%% crystallography ---------------------------------------------------------
% One crystalSymmetry object drives the lath crystallography, the orientation
% grid AND the IPF colouring, so every stage shares one convention.
% a = 2.95 A, c = 4.68 A -> c/a = 1.586.
exp_para.crystal = crystalSymmetry('6/mmm',[2.95 2.95 4.68], ...
    'X||a*','Y||b','mineral','Ti-alpha');

%% build the alpha-lath etching model -------------------------------------
% Every parameter below is a knob on the mechanism described at the top of this
% file; the defaults are the physically-motivated ones.  See lathGeometryTi64
% for the derivation of each, and lathReflectors for where it enters.
%
% WHICH CRYSTAL FEATURES FORM THE CORRUGATION
%   lathHabit    "rectangular" (default) treats the alpha lath as the 6-faced
%                platelet it actually is: 2 basal broad faces, 2 PRISMATIC
%                side faces and 2 side faces containing the LATTICE INVARIANT
%                LINE (the growth direction, Fig 2.14).  Its basal face is
%                then an approximate rectangle with just TWO rim directions.
%                "hexagonal" offers all six <10-10> equivalents instead, which
%                is an 8-faced hexagonal platelet - kept only for comparison.
%   lilUVW       lattice invariant line in the PARENT BETA frame; [-3 3 5]
%                from Fig 2.14b.  It is carried into the alpha frame across
%                the Burgers OR using that figure's own beta variant
%                ((110)beta // (0001)alpha), NOT variant 1 of Table 2.1 - the
%                setup checks that the result lands in the basal plane and
%                reports its angle to <11-20> (14.42 deg, matching the figure).
%   edgeWeights  [prismatic rim, lattice-invariant-line rim] relative weight -
%                effectively the lath ASPECT RATIO.  A lath grows long along
%                the LIL and narrow across it, so the two LIL rims are far
%                longer than the two prismatic end rims.  Use [1 1] for an
%                equiaxed platelet.  The two rim directions are 75.6 deg
%                apart, so one lath alone can give a NON-ORTHOGONAL pair of
%                stripes - an alternative to variant superposition.
%   rimBias      -1..1, relative weight of the two PARALLEL rims of the basal
%                rectangle (default 0 = equal).  The two are not physically
%                equivalent - beta sits above one and below the other - but
%                the chapter does not quantify it.  Non-zero makes the
%                parallel-to-c stripe one-sided.
%   ledgeUVTW    which prismatic side faces are allowed, via the basal-plane
%                direction they cut.  The chapter names [10-10] as the ledge
%                terrace direction at point 4.  Use [1 0 -1 0; 1 1 -2 0] to
%                allow <11-20> rims as well.
%   ledgeSelect  "best"  one realised edge per lath (variant selection),
%                        phi2-sensitive - the default;
%                "trace" the exact basal surface trace, so the stripe is
%                        exactly parallel to c and phi2 does nothing;
%                "all"   every equivalent -> 6-spoke asterisk (diagnostic).
%   prismHKIL    prism planes forming the ledge terraces; the pair named in
%                Fig 2.15 is {10-10} + {11-20}.
%   prismSelect  "all" lets the theta window select which terrace is visible.
%
% RELATIVE BRIGHTNESS OF THE FOUR FEATURES
%   ampTerraceRidge  stripe PERPENDICULAR to c  (point-1 mechanism)
%   ampBasalEdge     stripe PARALLEL to c       (point-4 mechanism)
%   ampBasalPlateau  basal plane-mirror point   (theta = 90 - 2*Phi)
%   ampPrismPlateau  prism ledge-terrace point  (the R2 of Fig 2.15)
%   pyramidFrac      stripe -> point trade when the terrace ridge pierces the
%                    surface instead of lying in it (the phi2 effect at the
%                    end of section 2.3.2)
%
% SHAPE / WEIGHTING
%   exposeExp    exponent on the cos(Phi) / sin(Phi) exposure factors
%   coincExp     exponent on cos(tilt): how strictly a crest line has to lie
%                in the surface.  Raise it to sharpen the phi2 sensitivity.
%   stripeWidth  cross-stripe half-width, deg
%   peakWidth    specular point half-width, deg
%   bevelSpan    EXTRA arc beyond the facet wedge, deg (default 0).  A real
%                crest is a bevel joining two facets, and rounding it sweeps
%                the normal from one face round to the other and no further -
%                so the lit arc IS the wedge angle, with no free parameter:
%                arcSpan = arcWedge + bevelSpan, centred on the bisector of
%                the two facets.  arcWedge is 30 deg for the prism-prism
%                terrace ridge and exactly 90 deg for a basal rim.  Raise
%                bevelSpan only if a measured stripe is longer than the wedge
%                allows; >= 360 switches the gate off (symmetric stripe).
%   bevelBias    signed tilt of the TERRACE RIDGE arc centre, deg.  The two
%                prism facets that crest joins straddle the surface normal, so
%                the model predicts a nearly symmetric perpendicular-to-c
%                stripe.  The real crest is asymmetric because one flank is
%                the etch-resistant beta lamella, but that flank is not
%                derivable from the alpha orientation alone - so fit this.
%                Expect order 10-20 deg; the sign picks which lobe brightens.
%   bevelSoft    soft shoulder of that arc, deg (0 = hard cut)
%   stripeTaper  exponent on (n.Z) along the stripe; 2 is the microfacet
%                foreshortening pair (illumination x visible area), which
%                brightens the stripe toward high theta as observed
%   shadowTh     specular points below this theta are damped by shadowing
%   combine      "max" brightest reflector wins | "sum" incoherent addition
%
% BURGERS-OR SIBLING VARIANTS (experimental, off by default)
%   variantSuperpose  true adds sibling-variant stripes, the chapter's
%                     explanation for non-orthogonal stripes in one DRM pixel
%   variantWeight     relative weight of those siblings
exp_para = lathGeometryTi64(exp_para, ...
    lathHabit       = "rectangular", ...
    lilUVW          = [-3 3 5], ...
    edgeWeights     = [1 1], ... % they are not fully dug out by etching
    rimBias         = 0, ...
    ledgeUVTW       = [1 0 -1 0], ...
    ledgeSelect     = "best", ...
    prismHKIL       = [1 0 -1 0; 1 1 -2 0], ...
    prismSelect     = "all", ...
    ampTerraceRidge = 1.00, ...
    ampBasalEdge    = 1.00, ...
    ampBasalPlateau = 0.25, ...
    ampPrismPlateau = 0.25, ...
    pyramidFrac     = 0.50, ...
    exposeExp       = 1, ...
    coincExp        = 2, ...
    stripeWidth     = 6, ...
    peakWidth       = 7, ...
    bevelSpan       = 0, ...
    bevelBias       = 0, ...
    bevelSoft       = 30, ...
    stripeTaper     = 2, ...
    shadowTh        = 20, ...
    combine         = "max", ...
    variantSuperpose = false, ...
    variantWeight    = 0.35, ...
    selfCheck        = true);

%% sanity check the mechanism with NO data --------------------------------
% Reproduces the four Fig 2.12 points before any file paths are involved:
% expect a stripe perpendicular to c to dominate for points 1-3 (high Phi) and
% a stripe parallel to c to dominate for point 4 (Phi = 37.8 deg), plus a
% basal specular point at theta = 90 - 2*Phi = 14.4 deg for point 4.
euPaper = [211.4 87.9  4.8;    % Fig 2.12 point 1
           160.3 87.7  4.2;    % Fig 2.12 point 2
           180.0 72.4  3.9;    % Fig 2.12 point 3
           289.2 37.8 28.0];   % Fig 2.12 point 4
checkTi64Mechanism(euPaper, exp_para, ...
    labels = ["Fig2.12 pt1";"Fig2.12 pt2";"Fig2.12 pt3";"Fig2.12 pt4"]);

% The stripes should NOT be symmetric about the DRP centre.  The printed
% brightAz/brightTh columns say where each bevel puts its maximum, and the
% circles on the overlay mark it.  To see how much of the asymmetry is the
% bevel, compare against the idealised full semi-cylinder:
%   epSym = lathGeometryTi64(exp_para, bevelSpan=360);
%   checkTi64Mechanism(euPaper, epSym, verbose=false);
%
% The two bevel knobs behave independently and are the ones to FIT against
% measured DRPs:
%   bevelSpan controls HOW asymmetric a stripe is.  On point 4 the ratio of
%     the bright lobe to the far lobe goes 20 deg -> 35:1, 60 -> 13:1,
%     90 -> 2:1, 150 -> 1:1 (symmetric).
%   bevelBias controls the perpendicular-to-c stripe specifically, which is
%     otherwise predicted near-symmetric because its two prism facets straddle
%     the surface normal.  Sweep it against a high-Phi grain:
%       for b = -20:10:20
%           ep = lathGeometryTi64(exp_para, bevelBias=b);
%           checkTi64Mechanism(euKnown, ep, measured=measDRP, verbose=false);
%       end

% How each Euler angle should behave (section 2.3.2, closing bullets):
%   phi1 -> pure azimuthal rotation of the whole DRP (c-axis azimuth = phi1-90)
%   Phi  -> swaps which stripe dominates, and sets the specular theta = 90-2Phi
%   phi2 -> weak: rotates the cell about c, moving the realised ledge edge off
%           the exact basal trace (dims the || stripe and swings its azimuth)
% Sweep them to confirm:
%   checkTi64Mechanism([0 20 0; 0 40 0; 0 60 0; 0 80 0], exp_para);
%   checkTi64Mechanism([0 60 0; 0 60 10; 0 60 20; 0 60 30], exp_para);

%% load sample dataset ----------------------------------------------------
% <<< PLUG IN YOUR PATHS HERE >>>
% drp_loader takes the image folder through its 'folder' option; the default
% inside drp_loader still points at the cpTi dataset.
%   [igray_sample, phitheta, pos, img_sample] = drp_loader(exp_para, pos1, ...
%       format='jpg', scale=scaleCoeff, folder="C:\path\to\Ti64\DRM\images");
[igray_sample, phitheta, pos, img_sample] = drp_loader( ...
    exp_para, pos1, format='jpg', scale=scaleCoeff);

% convert into a DRP stack: [n1 x n2] cells, each th_num x ph_num
drp_original = igray2drp(igray_sample, phitheta, exp_para);

%% verify the mechanism against MEASURED DRPs -----------------------------
% Pick a handful of grains whose EBSD orientation you know, read their DRPs out
% of drp_original, and compare.  NOTE the DRP viewer reports (col,row) while
% drp_original is indexed {row,col}.
%
% <<< PLUG IN YOUR EBSD ORIENTATIONS AND MATCHING PIXELS HERE >>>
euKnown = euPaper;                     % K x 3 [phi1 Phi phi2] from EBSD, deg
% m1 = drp_original{row1, col1};
% m2 = drp_original{row2, col2};
% m3 = drp_original{row3, col3};
% m4 = drp_original{row4, col4};
% measDRP = {m1; m2; m3; m4};
%
% checkTi64Mechanism(euKnown, exp_para, measured=measDRP);

% The DRM azimuth zero and the EBSD sample X are mounted at a fixed rotation
% about the sample normal, and the chapter notes a residual rotational
% misalignment.  Sweep phiOffset until ONE value aligns the predicted stripes
% for ALL grains - that value is your EBSD -> DRM registration.  It does not
% affect indexing (the engine absorbs rotation about Z), only this check.
% for off = 0:15:165
%     checkTi64Mechanism(euKnown, exp_para, measured=measDRP, ...
%         phiOffset=off, verbose=false);
% end

% Once registered, tune the model on the same grains, e.g.
%   ep = lathGeometryTi64(exp_para, coincExp=4, stripeWidth=4, ...
%                         ampBasalPlateau=0.1, ampPrismPlateau=0.1);
%   checkTi64Mechanism(euKnown, ep, measured=measDRP, phiOffset=<yours>);
% and, to test the chapter's variant-superposition explanation of
% non-orthogonal stripes:
%   ep = lathGeometryTi64(exp_para, variantSuperpose=true, variantWeight=0.4);
%   checkTi64Mechanism(euKnown, ep, measured=measDRP, phiOffset=<yours>);

%% DRP dictionary ---------------------------------------------------------
% CAVEAT: a lath DRP is two orthogonal stripes through the centre, so its
% azimuth signature has 180 deg period -> phi1 is recovered modulo 180.  phi2
% is only weakly encoded (through the ledge-edge selection), so expect the
% indexing to constrain the c-axis well and phi2 poorly.
if ~use_saved_drp_dic
    [drpDic, euDic, rotDic] = makeDRPdic_Ti64(exp_para, resolution=3);
    save('DRP_dictionary_Ti64_lath.mat','drpDic','euDic','rotDic');
else
    S = load('DRP_dictionary_Ti64_lath.mat','drpDic','euDic','rotDic');
    drpDic = S.drpDic; euDic = S.euDic; rotDic = S.rotDic;
end

%% autoencoder ------------------------------------------------------------
if ~use_saved_ae
    AE_DRM = trainAutoencoder(drpDic, 100, ...
        'MaxEpochs',200, 'L2WeightRegularization',0.001, ...
        'SparsityRegularization',4, 'SparsityProportion',0.10, ...
        'ScaleData',false, 'UseGPU',false);
    save('DRP_autoencoder_Ti64_lath.mat','AE_DRM');
else
    S = load('DRP_autoencoder_Ti64_lath.mat','AE_DRM');
    AE_DRM = S.AE_DRM;
end

%% index ------------------------------------------------------------------
tic
index_result = IndexingEngine(drp_original, AE_DRM, exp_para, drpDic, euDic, rotDic);
toc

%% results ---------------------------------------------------------------
[n1,n2]   = size(drp_original);
index_num = cellfun(@(x) sum(x,'all'), drp_original);
non_index_bg = index_num > 3e4;                    % tune this threshold
figure, imshow(non_index_bg); title('background mask (tune threshold)');

figure, imshow(plot_ipf_map(index_result.EUmap));
title('IPF-z, wrought Ti64 alpha (lath mechanism)');

%% validate: measured vs predicted DRP at clicked pixels -----------------
[drp_measurement, drp_predicted, xy] = check_indexing_result_Ti64( ...
    index_result.EUmap, drp_original, exp_para);

%% (optional) direct indexing without the autoencoder --------------------
% drpLib     = DRPLibGenerator_Ti64(5*degree, exp_para);
% indexResult = DirectDIEngine(drp_original, drpLib);
% figure, imshow(plot_ipf_map(indexResult.euMap));
% title('direct indexing (no autoencoder), Ti64 alpha lath mechanism');
