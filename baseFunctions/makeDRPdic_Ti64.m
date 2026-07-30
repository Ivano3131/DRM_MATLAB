function [drpDic, euDic, rotDic] = makeDRPdic_Ti64(exp_para, options)
% makeDRPdic_Ti64  Simulated DRP dictionary for wrought-Ti64 alpha indexing.
%
%   [drpDic, euDic, rotDic] = makeDRPdic_Ti64(exp_para)
%   [...] = makeDRPdic_Ti64(exp_para, resolution=3, verbose=true)
%
% Identical bookkeeping to makeDRPdic_hcp, but the patterns come from
% DRPsim_Ti64 (the alpha-lath etching mechanism) instead of DRPsim_hcp.
%
% Samples the hexagonal fundamental zone with equispacedSO3Grid (using
% exp_para.crystal, so the symmetry matches the IPF colouring), simulates a
% DRP per orientation, and circularly shifts each one in azimuth so its
% strongest column sits at phi = 0.  The applied shift is returned in rotDic
% and folded back into phi1 during indexing, which makes matching invariant
% to rotation of the grain about the sample Z axis.
%
% CAVEAT specific to this mechanism --------------------------------------
% A lath DRP is essentially two orthogonal stripes through the centre, so its
% column-sum profile has 180 deg period: the azimuth alignment above is
% 180 deg ambiguous and phi1 comes out modulo 180 rather than 360.  Expect
% that in the indexed map, and note that IndexingEngine aligns measured DRPs
% with a different statistic (count of maximum-valued pixels per column) than
% the one used here (column sum) - worth harmonising before trusting phi1.
%
% Outputs
%   drpDic  n x 1 cell, each a th_num x ph_num double DRP in [0,1]
%   euDic   n x 3   Bunge Euler angles [phi1 Phi phi2] in degrees
%   rotDic  n x 1   azimuth shift (in ph_num steps) applied to each DRP
% -------------------------------------------------------------------------
arguments
    exp_para struct
    options.resolution (1,1) double = 3      % orientation grid spacing, deg
    options.verbose (1,1) logical = true
end

if ~isfield(exp_para,'lath') || isempty(exp_para.lath)
    error('makeDRPdic_Ti64:notSetup', ...
        'exp_para.lath is missing. Run exp_para = lathGeometryTi64(exp_para) first.');
end

cs  = exp_para.crystal;
ori = equispacedSO3Grid(cs,'resolution',options.resolution*degree);
nn  = length(ori);

% pull Euler angles out as plain arrays so the parfor body needs no MTEX
phi1 = ori.phi1 ./ degree;
Phi  = ori.Phi  ./ degree;
phi2 = ori.phi2 ./ degree;

drpDic = cell(nn,1);
euDic  = zeros(nn,3);
rotDic = zeros(nn,1);

% strip the MTEX object; exp_para.lath is pure numeric/string data
ep = exp_para;
ep.crystal = [];

parfor index = 1:nn
    tmpDRP = DRPsim_Ti64(phi1(index), Phi(index), phi2(index), ep);
    colsum = sum(tmpDRP,1);
    [~,shift] = max(colsum);
    shift = -shift;
    drpDic{index}  = double(circshift(tmpDRP,shift,2)) / 255;
    euDic(index,:) = [mod(phi1(index),360), Phi(index), mod(phi2(index),360)];
    rotDic(index)  = shift;
end

if options.verbose
    fprintf('Ti64 lath DRP dictionary built: %d orientations at %.1f deg resolution.\n', ...
        nn, options.resolution);
end
end
