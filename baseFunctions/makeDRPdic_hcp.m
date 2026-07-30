function [drpDic, euDic, rotDic] = makeDRPdic_hcp(exp_para, options)
% makeDRPdic_hcp  Build the simulated DRP dictionary for HCP indexing.
%
%   [drpDic, euDic, rotDic] = makeDRPdic_hcp(exp_para)
%   [...] = makeDRPdic_hcp(exp_para, resolution=3, verbose=true)
%
% Samples the hexagonal orientation fundamental zone with equispacedSO3Grid
% (using exp_para.crystal, so the symmetry matches the IPF colouring), then
% simulates a DRP for each orientation with DRPsim_hcp.  Each DRP is circularly
% shifted in azimuth so its dominant column sits at phi = 0; the applied shift
% is stored in rotDic and the residual azimuth is folded back during indexing,
% making the match invariant to rotation of the grain about the sample Z axis.
%
% Outputs
%   drpDic  goodnn x 1 cell, each a th_num x ph_num double DRP in [0,1]
%   euDic   goodnn x 3   Bunge Euler angles [phi1 Phi phi2] in degrees
%   rotDic  goodnn x 1   azimuth shift (in ph_num steps) applied to each DRP
% -------------------------------------------------------------------------
arguments
    exp_para struct
    options.resolution (1,1) double = 3      % orientation grid spacing, deg
    options.verbose (1,1) logical = true
end

if ~isfield(exp_para,'facetNormals')
    exp_para = facetNormalsHCP(exp_para);
end
cs = exp_para.crystal;

ori = equispacedSO3Grid(cs,'resolution',options.resolution*degree);
goodnn = length(ori);

% pull Euler angles out as plain arrays so the parfor body needs no MTEX object
phi1 = ori.phi1 ./ degree;
Phi  = ori.Phi  ./ degree;
phi2 = ori.phi2 ./ degree;

drpDic = cell(goodnn,1);
euDic  = zeros(goodnn,3);
rotDic = zeros(goodnn,1);

% strip the MTEX object so each worker carries only numeric data
ep = exp_para;
ep.crystal = [];

parfor index = 1:goodnn
    tmpDRP = DRPsim_hcp(phi1(index),Phi(index),phi2(index),ep);
    colsum = sum(tmpDRP,1);
    [~,shift] = max(colsum);
    shift = -shift;
    drpDic{index} = double(circshift(tmpDRP,shift,2)) / 255;
    euDic(index,:) = [mod(phi1(index),360), Phi(index), mod(phi2(index),360)];
    rotDic(index) = shift;
end

if options.verbose
    fprintf('HCP DRP dictionary built: %d orientations at %.1f deg resolution.\n', ...
        goodnn, options.resolution);
end
end
