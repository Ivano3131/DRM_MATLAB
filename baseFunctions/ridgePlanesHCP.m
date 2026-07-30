function exp_para = ridgePlanesHCP(exp_para)
% ridgePlanesHCP  Build crystallographic PLANE normals whose SURFACE TRACES
%                 form the ridges (for the ridge/stripe DRP model).
%
%   exp_para = ridgePlanesHCP(exp_para)
%
% Physical model: the smooth elongated ridges seen by AFM/SEM are the surface
% traces of a crystallographic plane.  A plane with normal m intersects the
% sample surface (normal Z) along the line t = m x Z, which lies IN the surface
% (horizontal) for ANY grain orientation.  A cylindrical ridge along t reflects
% the fixed camera along the great circle perpendicular to t; because t is
% horizontal that great circle passes through +Z, so the DRP stripe ALWAYS goes
% through the pattern centre (as observed).  The stripe azimuth is set by the
% in-plane direction of m, i.e. by the crystallography + grain orientation.
%
% The trace t = m x Z depends on orientation, so it is formed per-orientation
% inside DRPsim_hcp; this function only builds the crystal-frame plane normals m
% (symmetry-equivalent set, correct hexagonal metric) that get rotated there.
%
% Required input field
%   exp_para.ridgeHKIL   N x 4 matrix, each row a Miller-Bravais plane {h k i l}
%                        (i = -(h+k)) whose traces form the ridges, e.g.
%                          [0 0  0 1]  basal (0001)       -> one trace / stripe
%                          [1 0 -1 0]  prism {10-10}
%                          [1 0 -1 1]  pyramidal {10-11}
%
% Optional input fields
%   exp_para.crystal        MTEX crystalSymmetry (default: alpha-Ti 6/mmm)
%   exp_para.ridgeWeights   N x 1 relative weight per plane family
%
% Output fields added to exp_para
%   exp_para.ridgePlaneNormals  M x 3 unit plane normals in the crystal frame
%   exp_para.ridgeW             M x 1 weights (sum-normalised)
%   exp_para.ridgeFamilyId      M x 1 family index (1..N) per normal
% -------------------------------------------------------------------------
arguments
    exp_para struct
end

if ~isfield(exp_para,'crystal') || isempty(exp_para.crystal)
    exp_para.crystal = crystalSymmetry('6/mmm',[2.95 2.95 4.68], ...
        'X||a*','Y||b','mineral','Ti-alpha');
end
cs = exp_para.crystal;

if ~isfield(exp_para,'ridgeHKIL') || isempty(exp_para.ridgeHKIL)
    error('ridgePlanesHCP:noFamily', ...
        'exp_para.ridgeHKIL is required (N x 4 list of {h k i l} ridge planes).');
end
hkil = exp_para.ridgeHKIL;
if size(hkil,2) ~= 4
    error('ridgePlanesHCP:badFamily','ridgeHKIL must have 4 columns (h k i l).');
end
N = size(hkil,1);

if isfield(exp_para,'ridgeWeights') && ~isempty(exp_para.ridgeWeights)
    rw = exp_para.ridgeWeights(:);
else
    rw = ones(N,1);
end

normals = zeros(0,3);
ridgeW  = zeros(0,1);
famId   = zeros(0,1);

for k = 1:N
    m  = Miller(hkil(k,1),hkil(k,2),hkil(k,3),hkil(k,4),cs,'hkl');
    ms = symmetrise(m);
    xyz = normr(ms.xyz);
    xyz = local_canonicalSign(xyz);
    xyz = uniquetol(xyz,1e-6,'ByRows',true);

    nk = size(xyz,1);
    normals = [normals; xyz];             %#ok<AGROW>
    ridgeW  = [ridgeW;  rw(k)*ones(nk,1)]; %#ok<AGROW>
    famId   = [famId;   k*ones(nk,1)];    %#ok<AGROW>
end

if sum(ridgeW) > 0, ridgeW = ridgeW / sum(ridgeW); end

exp_para.ridgePlaneNormals = normals;
exp_para.ridgeW            = ridgeW;
exp_para.ridgeFamilyId     = famId;
end

% -------------------------------------------------------------------------
function xyz = local_canonicalSign(xyz)
tol = 1e-8;
for r = 1:size(xyz,1)
    v = xyz(r,:);
    idx = find(abs(v) > tol,1,'first');
    if ~isempty(idx) && v(idx) < 0
        xyz(r,:) = -v;
    end
end
end
