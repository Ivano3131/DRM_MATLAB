function out = cubeAmbiguity(exp_para, options)
% cubeAmbiguity  Can this cube-lath model be indexed at all?
%
%   out = cubeAmbiguity(exp_para)
%   out = cubeAmbiguity(exp_para, resolution=5, nTest=60, corrTol=0.995)
%
% ASK THIS BEFORE TRUSTING A MAP.  A DRP of stripes through the pattern centre
% carries only so much information, and some box shapes throw information away
% outright: if two of the three axes have the same length and the same edge
% amplitude, a 90 deg rotation about the third maps the box onto itself, and
% two genuinely different crystal orientations produce the SAME DRP.  A
% literal cube (dims [1 1 1], edgeAmp [1 1 1]) has the full octahedral
% symmetry and is completely unindexable.
%
% This measures both halves of the question:
%
%   AMBIGUITY RADIUS - for each dictionary entry, the largest misorientation
%     to any OTHER entry (over all azimuth shifts) that still correlates above
%     corrTol.  Near 0 means that orientation is pinned down; 90 or 120 deg
%     means it has a far-away twin the matcher cannot rule out.  The p90 is
%     the number to watch: it says what fraction of orientation space is
%     hopeless.
%
%   ROUND TRIP - simulate a random orientation, match it back, and report the
%     misorientation.  The median should be BELOW the grid step; the tail is
%     the same degeneracy showing up in practice.
%
% Both are measured modulo the box's own proper symmetry (identity plus 180
% deg about each box axis), which is the group the model genuinely cannot and
% should not resolve.
%
% OPTIONS -----------------------------------------------------------------
%   resolution  dictionary spacing, deg (default 5)
%   nTest       random orientations for the round trip (default 60)
%   corrTol     correlation above which two patterns count as indistinguishable
%               (default 0.995)
%   verbose     true (default) prints the summary
%
% OUTPUT ------------------------------------------------------------------
%   out.radius     n x 1 ambiguity radius per dictionary entry, deg
%   out.roundTrip  nTest x 1 round-trip misorientation, deg
%   out.euDic      the dictionary Euler pairs, so a bad region can be located
%   out.verdict    string
%
% See also cubeGeometryTi64, makeDRPdic_cube, matchDRPcube.
% -------------------------------------------------------------------------
arguments
    exp_para struct
    options.resolution (1,1) double = 5
    options.nTest (1,1) double = 60
    options.corrTol (1,1) double = 0.995
    options.verbose (1,1) logical = true
end

if ~isfield(exp_para,'cube') || isempty(exp_para.cube)
    error('cubeAmbiguity:notSetup', ...
        'exp_para.cube is missing. Run exp_para = cubeGeometryTi64(exp_para) first.');
end

[drpDic, euDic] = makeDRPdic_cube(exp_para, resolution=options.resolution, ...
    verbose=options.verbose);
nD     = numel(drpDic);
th_num = exp_para.th_num;
ph_num = exp_para.ph_num;

S = local_boxSym(exp_para.cube);

% ---- correlate the dictionary against itself, over all shifts ----------
D = zeros(th_num, ph_num, nD);
for k = 1:nD
    d = drpDic{k};  d = d - mean(d(:));
    D(:,:,k) = d / max(norm(d,'fro'), eps);
end
A  = permute(conj(fft(D,[],2)), [3 1 2]);
Fm = permute(fft(D,[],2), [1 3 2]);
CC = real(ifft(pagemtimes(A,Fm), [], 3));       % nD x nD x shift

radius = zeros(nD,1);
for k = 1:nD
    cc = squeeze(CC(:,k,:));
    Rk = local_bunge(0, euDic(k,1), euDic(k,2));
    [hit, sh] = find(cc > options.corrTol);
    a = 0;
    for h = 1:numel(hit)
        Rj = local_bunge((sh(h)-1)*(360/ph_num), euDic(hit(h),1), euDic(hit(h),2));
        a  = max(a, local_misAngle(Rk, Rj, S));
    end
    radius(k) = a;
end

% ---- round trip --------------------------------------------------------
rt = zeros(options.nTest,1);
for ii = 1:options.nTest
    t = [360*rand, 90*rand, 180*rand];
    r = matchDRPcube(DRPsim_cube(t(1),t(2),t(3),exp_para), drpDic, euDic, ...
        exp_para, verbose=false);
    rt(ii) = local_misAngle(local_bunge(t(1),t(2),t(3)), ...
                            local_bunge(r.euler(1),r.euler(2),r.euler(3)), S);
end

p90 = prctile(radius,90);
if p90 < 2*options.resolution
    verdict = "INDEXABLE - the dictionary is essentially unique.";
elseif p90 < 30
    verdict = "USABLE - a minority of orientations have close twins.";
else
    verdict = "DEGENERATE - a large part of orientation space has far-away " + ...
              "twins with the same DRP. Give the three box axes distinct " + ...
              "dims (or distinct edgeAmp) before trusting a map.";
end

out.radius    = radius;
out.roundTrip = rt;
out.euDic     = euDic;
out.verdict   = verdict;

if options.verbose
    fprintf(['\ncubeAmbiguity  (dims [%.2f %.2f %.2f], etchFrac %.2f, ', ...
             '%d entries at %g deg)\n'], exp_para.cube.dims, ...
             exp_para.cube.etchFrac, nD, options.resolution);
    fprintf('  ambiguity radius (deg): median %5.1f   p90 %5.1f   max %5.1f\n', ...
        median(radius), p90, max(radius));
    fprintf('  round trip       (deg): median %5.1f   p90 %5.1f   |  %d/%d beyond 10 deg\n', ...
        median(rt), prctile(rt,90), sum(rt>10), numel(rt));
    fprintf('  -> %s\n', verdict);
end
end

% =========================================================================
function S = local_boxSym(cube)
% The box's own proper symmetry: identity plus 180 deg about each box axis
% (point group 222).  This is what the model legitimately cannot resolve, so
% misorientations are measured modulo it.
B = [cube.eC; cube.eB; cube.eP];
S = cell(1,4);
S{1} = eye(3);
for k = 1:3
    a = B(k,:);
    S{k+1} = 2*(a.'*a) - eye(3);
end
end

% =========================================================================
function R = local_bunge(p1, P, p2)
% Same convention as EulerRotate: v_sample = v_crystal * R.
Rz = @(t)[cosd(t) sind(t) 0; -sind(t) cosd(t) 0; 0 0 1];
Rx = @(t)[1 0 0; 0 cosd(t) sind(t); 0 -sind(t) cosd(t)];
R  = Rz(p2) * Rx(P) * Rz(p1);
end

% =========================================================================
function a = local_misAngle(R1, R2, S)
a = 180;
for ii = 1:numel(S)
    M = R1 * (S{ii}*R2).';
    a = min(a, acosd(min(max((trace(M)-1)/2,-1),1)));
end
end
