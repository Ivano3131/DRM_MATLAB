function exp_para = facetNormalsHCP(exp_para)
% facetNormalsHCP  Build crystallographically-correct HCP facet normals.
%
%   exp_para = facetNormalsHCP(exp_para)
%
% Turns a list of {hkil} plane families into the full set of symmetry-
% equivalent facet-plane normals, expressed in the CRYSTAL Cartesian frame,
% using MTEX so that the hexagonal metric (a, c -> correct c/a ratio) and the
% 6/mmm symmetry are applied correctly.  These normals are consumed by
% rotate_facet_hcp / DRPsim_hcp and are consistent with the orientation grid
% (equispacedSO3Grid) and the IPF colouring (ipfHSVKey), because all three use
% the SAME crystalSymmetry object.
%
% Required input fields
%   exp_para.facetHKIL   N x 4 matrix, each row a Miller-Bravais plane {h k i l}
%                        with i = -(h+k).  Examples:
%                          [1 0 -1 0]  first-order prism {10-10}
%                          [1 1 -2 0]  second-order prism {11-20}
%                          [1 0 -1 1]  first-order pyramidal {10-11}
%                          [0 0  0 1]  basal (0001)
%                        Multiple rows = multiple co-existing facet families.
%
% Optional input fields
%   exp_para.crystal          MTEX crystalSymmetry (default: alpha-Ti 6/mmm)
%   exp_para.facetFaceWeights N x 1 relative weight of each family's PEAKS
%   exp_para.facetBandWeights N x 1 relative weight of each family's BANDS
%
% Output fields added to exp_para
%   exp_para.facetNormals    M x 3 unit normals in the crystal Cartesian frame
%   exp_para.faceW           M x 1 peak weights (sum-normalised)
%   exp_para.pairW           M x 1 band weights (sum-normalised)
%   exp_para.facetFamilyId   M x 1 index (1..N) of the family each normal came from
%
% -------------------------------------------------------------------------
arguments
    exp_para struct
end

% default alpha-Ti crystal symmetry (a = 2.95 A, c = 4.68 A -> c/a = 1.586)
if ~isfield(exp_para,'crystal') || isempty(exp_para.crystal)
    exp_para.crystal = crystalSymmetry('6/mmm',[2.95 2.95 4.68], ...
        'X||a*','Y||b','mineral','Ti-alpha');
end
cs = exp_para.crystal;

if ~isfield(exp_para,'facetHKIL') || isempty(exp_para.facetHKIL)
    error('facetNormalsHCP:noFamily', ...
        'exp_para.facetHKIL is required (N x 4 list of {h k i l} plane families).');
end
hkil = exp_para.facetHKIL;
if size(hkil,2) ~= 4
    error('facetNormalsHCP:badFamily','facetHKIL must have 4 columns (h k i l).');
end
N = size(hkil,1);

% per-family relative weights (default uniform)
if isfield(exp_para,'facetFaceWeights') && ~isempty(exp_para.facetFaceWeights)
    fw = exp_para.facetFaceWeights(:);
else
    fw = ones(N,1);
end
if isfield(exp_para,'facetBandWeights') && ~isempty(exp_para.facetBandWeights)
    bw = exp_para.facetBandWeights(:);
else
    bw = ones(N,1);
end

normals = zeros(0,3);
faceW   = zeros(0,1);
pairW   = zeros(0,1);
famId   = zeros(0,1);

for k = 1:N
    % define the plane and generate all symmetry-equivalent normals.
    % 6/mmm is centrosymmetric, so symmetrise already includes -n; the
    % canonical-sign + uniquetol step below collapses each +/- pair to one
    % distinct facet plane.
    m  = Miller(hkil(k,1),hkil(k,2),hkil(k,3),hkil(k,4),cs,'hkl');
    ms = symmetrise(m);                 % all symmetry-equivalent normals
    xyz = normr(ms.xyz);                % Cartesian plane normals, unit length
    xyz = local_canonicalSign(xyz);
    xyz = uniquetol(xyz,1e-6,'ByRows',true);

    nk = size(xyz,1);
    normals = [normals; xyz];                %#ok<AGROW>
    faceW   = [faceW;  fw(k)*ones(nk,1)];    %#ok<AGROW>
    pairW   = [pairW;  bw(k)*ones(nk,1)];    %#ok<AGROW>
    famId   = [famId;  k*ones(nk,1)];        %#ok<AGROW>
end

if sum(faceW) > 0, faceW = faceW / sum(faceW); end
if sum(pairW) > 0, pairW = pairW / sum(pairW); end

exp_para.facetNormals   = normals;
exp_para.faceW          = faceW;
exp_para.pairW          = pairW;
exp_para.facetFamilyId  = famId;
end

% -------------------------------------------------------------------------
function xyz = local_canonicalSign(xyz)
% Flip each row so its first significant component is positive, so that a
% normal and its antipode map to the same representative row.
tol = 1e-8;
for r = 1:size(xyz,1)
    v = xyz(r,:);
    idx = find(abs(v) > tol,1,'first');
    if ~isempty(idx) && v(idx) < 0
        xyz(r,:) = -v;
    end
end
end
