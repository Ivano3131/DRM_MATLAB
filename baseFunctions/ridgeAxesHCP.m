function exp_para = ridgeAxesHCP(exp_para)
% ridgeAxesHCP  Build crystallographic ridge (cylinder) axes for the HCP
%               "ridge reflector" DRP model.
%
%   exp_para = ridgeAxesHCP(exp_para)
%
% Companion to facetNormalsHCP, but for the smooth ELONGATED RIDGE topology
% seen by AFM/SEM (half-pipe / cylindrical features) rather than flat facets.
% A cylindrical ridge with axis t has surface normals filling the great circle
% perpendicular to t, so it reflects the fixed overhead camera along that whole
% great circle -> a STRIPE in the DRP whose pole is the ridge axis (see
% DRPsim_hcp, reflectorModel = "ridge").
%
% This turns a list of crystallographic ridge DIRECTIONS <uvtw> into the full
% symmetry-equivalent set of unit axes in the CRYSTAL Cartesian frame, using
% MTEX so the hexagonal metric is applied correctly and everything stays on the
% same convention as the orientation grid and IPF map.
%
% Required input field
%   exp_para.ridgeUVTW   N x 4 matrix, each row a Miller-Bravais DIRECTION
%                        <u v t w> with u+v+t = 0.  Examples:
%                          [ 1  1 -2  0 ]  <11-20>  (a close-packed a-direction)
%                          [ 2 -1 -1  0 ]  <2-1-10> (an a-axis)
%                          [ 1  0 -1  0 ]  <10-10>
%                          [ 0  0  0  1 ]  <0001>   (c-axis)
%                          [ 1  1 -2  3 ]  <11-23>  (a pyramidal-ish direction)
%                        Multiple rows = multiple co-existing ridge orientations.
%
% Optional input fields
%   exp_para.crystal          MTEX crystalSymmetry (default: alpha-Ti 6/mmm)
%   exp_para.ridgeWeights     N x 1 relative weight of each direction family
%
% Output fields added to exp_para
%   exp_para.ridgeAxes        M x 3 unit ridge axes in the crystal Cartesian frame
%   exp_para.ridgeW           M x 1 weights (sum-normalised)
%   exp_para.ridgeFamilyId    M x 1 index (1..N) of the family each axis came from
% -------------------------------------------------------------------------
arguments
    exp_para struct
end

if ~isfield(exp_para,'crystal') || isempty(exp_para.crystal)
    exp_para.crystal = crystalSymmetry('6/mmm',[2.95 2.95 4.68], ...
        'X||a*','Y||b','mineral','Ti-alpha');
end
cs = exp_para.crystal;

if ~isfield(exp_para,'ridgeUVTW') || isempty(exp_para.ridgeUVTW)
    error('ridgeAxesHCP:noFamily', ...
        'exp_para.ridgeUVTW is required (N x 4 list of <u v t w> ridge directions).');
end
uvtw = exp_para.ridgeUVTW;
if size(uvtw,2) ~= 4
    error('ridgeAxesHCP:badFamily','ridgeUVTW must have 4 columns (u v t w).');
end
N = size(uvtw,1);

if isfield(exp_para,'ridgeWeights') && ~isempty(exp_para.ridgeWeights)
    rw = exp_para.ridgeWeights(:);
else
    rw = ones(N,1);
end

axes3   = zeros(0,3);
ridgeW  = zeros(0,1);
famId   = zeros(0,1);

for k = 1:N
    % 'uvw' -> lattice DIRECTION (uses the direct metric; c/a matters for
    % anything with a non-zero w).  A ridge axis is a line, so +/- are equal.
    m  = Miller(uvtw(k,1),uvtw(k,2),uvtw(k,3),uvtw(k,4),cs,'uvw');
    ms = symmetrise(m);
    xyz = normr(ms.xyz);
    xyz = local_canonicalSign(xyz);
    xyz = uniquetol(xyz,1e-6,'ByRows',true);

    nk = size(xyz,1);
    axes3  = [axes3;  xyz];              %#ok<AGROW>
    ridgeW = [ridgeW; rw(k)*ones(nk,1)]; %#ok<AGROW>
    famId  = [famId;  k*ones(nk,1)];    %#ok<AGROW>
end

if sum(ridgeW) > 0, ridgeW = ridgeW / sum(ridgeW); end

exp_para.ridgeAxes     = axes3;
exp_para.ridgeW        = ridgeW;
exp_para.ridgeFamilyId = famId;
end

% -------------------------------------------------------------------------
function xyz = local_canonicalSign(xyz)
% Flip each row so its first significant component is positive, so a direction
% and its antipode collapse to the same representative axis.
tol = 1e-8;
for r = 1:size(xyz,1)
    v = xyz(r,:);
    idx = find(abs(v) > tol,1,'first');
    if ~isempty(idx) && v(idx) < 0
        xyz(r,:) = -v;
    end
end
end
