function [drpDic, euDic] = makeDRPdic_cube(exp_para, options)
% makeDRPdic_cube  Simulated DRP dictionary for the Ti64 cube-lath model.
%
%   [drpDic, euDic] = makeDRPdic_cube(exp_para)
%   [drpDic, euDic] = makeDRPdic_cube(exp_para, resolution=5)
%
% WHY phi1 IS NOT SAMPLED -------------------------------------------------
% phi1 is a rotation about the sample Z axis, which is EXACTLY a rigid
% rotation of the DRP in azimuth - i.e. a circular shift of its columns.
% matchDRPcube searches all ph_num shifts anyway (for free, by FFT), so
% sampling phi1 would only duplicate patterns.  The dictionary is therefore
% built at phi1 = 0 and spans (Phi, phi2) only.
%
% WHY phi2 RUNS TO 180 AND NOT 60 ----------------------------------------
% The crystal is 6/mmm, but the BOX is not: a 60 deg rotation about c maps the
% {10-10} broad-face normal onto a different {10-10}, i.e. onto a DIFFERENT
% box.  The box's proper point group is 222, whose Euler fundamental zone is
%       phi1 in [0,360),  Phi in [0,90],  phi2 in [0,180).
% That is physically right: a real lath selects one variant, it does not
% average over the 6-fold.  The consequence is that the recovered phi2 is
% defined modulo 180 deg with respect to the box, and plot_ipf_map folds it
% back under 6/mmm for colouring.
%
% At the default 5 deg that is 19 x 36 = 684 patterns - small enough to build
% in seconds and to match by brute force, so there is no autoencoder and no
% kNN in this pipeline.  There is also no azimuth canonicalisation, which
% removes the dictionary-vs-engine shift-statistic mismatch documented in
% makeDRPdic_Ti64.m:20-23.
%
% INPUTS ------------------------------------------------------------------
%   exp_para  needs the DRP grid and .cube (cubeGeometryTi64).
%
% OPTIONS -----------------------------------------------------------------
%   resolution  grid spacing in DEGREES (default 5).
%   verbose     true (default) prints a summary.
%
% OUTPUTS -----------------------------------------------------------------
%   drpDic  n x 1 cell, each th_num x ph_num double in [0,1]
%   euDic   n x 2 [Phi phi2] in degrees (phi1 is 0 for every entry)
%
% See also DRPsim_cube, matchDRPcube, cubeGeometryTi64.
% -------------------------------------------------------------------------
arguments
    exp_para struct
    options.resolution (1,1) double {mustBePositive} = 5
    options.verbose (1,1) logical = true
end

if ~isfield(exp_para,'cube') || isempty(exp_para.cube)
    error('makeDRPdic_cube:notSetup', ...
        'exp_para.cube is missing. Run exp_para = cubeGeometryTi64(exp_para) first.');
end

res     = options.resolution;
PhiList = 0 : res : 90;
Ph2List = 0 : res : (180 - res);

[P2, P1] = meshgrid(Ph2List, PhiList);
Phi  = P1(:);
phi2 = P2(:);
n    = numel(Phi);

% MTEX objects cannot cross into parfor; nothing below needs the crystal.
ep = exp_para;
ep.crystal = [];

drpDic = cell(n,1);
euDic  = [Phi phi2];

if options.verbose
    fprintf('makeDRPdic_cube: %d patterns (Phi 0:%g:90, phi2 0:%g:%g) ...\n', ...
        n, res, res, 180-res);
end
tStart = tic;

parfor ii = 1:n
    drpDic{ii} = double(DRPsim_cube(0, Phi(ii), phi2(ii), ep)) / 255;
end

if options.verbose
    fprintf('makeDRPdic_cube: done in %.1f s.\n', toc(tStart));
end
end
