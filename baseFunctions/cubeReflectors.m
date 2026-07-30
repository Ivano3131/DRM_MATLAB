function E = cubeReflectors(eu1, eu2, eu3, cube)
% cubeReflectors  Rotate the lath box, cut it at the etch depth, list its edges.
%
%   E = cubeReflectors(eu1, eu2, eu3, cube)
%
% THE PHYSICS, fully exposed.  Every intermediate quantity is returned rather
% than folded into an amplitude, so the model can be audited number by number
% against a measured DRP.  Contains NO MTEX calls, so it is parfor-safe.
%
% WHAT IT DOES ------------------------------------------------------------
%   1. rotates the 8 box vertices and the 12 edge face-normals from the
%      crystal frame into the sample frame with EulerRotate (active Bunge
%      ZXZ: v_sample = v_crystal * Rz(phi1)*Rx(Phi)*Rz(phi2));
%   2. places the etch cut plane at
%           zcut = max(z) - etchFrac * (max(z) - min(z))
%      i.e. the top etchFrac of the ROTATED box height stands proud;
%   3. gives each edge an exposure = the fraction of its length above zcut.
%      z varies linearly along an edge, so this is exact, and it varies
%      smoothly with orientation - no on/off flicker between grid steps.
%
% INPUTS ------------------------------------------------------------------
%   eu1,eu2,eu3  Bunge Euler angles [phi1 Phi phi2], degrees.
%   cube         exp_para.cube from cubeGeometryTi64.
%
% OUTPUT ------------------------------------------------------------------
%   E  12x1 struct array, one per edge:
%     .family       "broadEnd" (axis || c) | "basalEnd" | "basalBroad"
%     .axisId       1 = || e_c, 2 = || e_b, 3 = || e_p
%     .axis         1x3 unit edge direction, SAMPLE frame, z-component >= 0
%     .n1 .n2       1x3 outward normals of the two faces the bevel joins
%     .bisector     1x3 normalise(n1+n2): the direction the bevel faces
%     .expose       0..1 fraction of the edge above the etch cut
%     .amp          edgeAmp * lengthWeight * expose^exposeExp
%     .tilt         deg, how far the edge tilts out of the surface plane.
%                   0 means the crest is horizontal, so its stripe is an exact
%                   DIAMETER through the DRP centre.
%     .stripeAz     deg mod 180, azimuth of that stripe
%     .arcCentreAz  deg 0..360, azimuth where the bevel faces
%     .arcCentreTh  deg, DRP theta of the same, = 2*(elevation(bisector) - 45).
%                   This is where the bevel puts its MAXIMUM along the stripe,
%                   i.e. the asymmetry a real rounded edge introduces.
%     .zMid         height of the edge midpoint (diagnostic)
%
% Also carries the per-orientation summary in every element (cheap, and it
% saves the caller recomputing it):
%     .cAxis        1x3 c-axis in the sample frame, z-component >= 0
%     .cAz          deg, c-axis azimuth (expect phi1 - 90)
%     .Phi          deg, c-axis polar angle folded to [0,90]
%
% See also cubeGeometryTi64, DRPsim_cube, checkCubeMechanism.
% -------------------------------------------------------------------------

nE = size(cube.edgeV,1);

% ---- rotate into the sample frame --------------------------------------
V  = EulerRotate(cube.vert,  eu1, eu2, eu3);      % 8x3
N1 = EulerRotate(cube.edgeN1, eu1, eu2, eu3);     % 12x3
N2 = EulerRotate(cube.edgeN2, eu1, eu2, eu3);     % 12x3

cS = EulerRotate(cube.eC, eu1, eu2, eu3);
if cS(3) < 0, cS = -cS; end                       % c is an axis, not a vector
cAz = mod(atan2d(cS(2), cS(1)), 360);
Phi = acosd(min(max(cS(3),-1),1));

% ---- the etch cut plane -------------------------------------------------
zTop = max(V(:,3));
zBot = min(V(:,3));
h    = zTop - zBot;
if h <= 0
    zcut = zTop - eps;            % degenerate: box lying exactly flat
else
    zcut = zTop - cube.etchFrac * h;
end

E = repmat(local_blank(), nE, 1);

for ii = 1:nE
    P1 = V(cube.edgeV(ii,1),:);
    P2 = V(cube.edgeV(ii,2),:);

    d = P2 - P1;
    L = norm(d);
    if L <= 0, continue; end
    t = d / L;
    if t(3) < 0, t = -t; end       % crest is an axis; fix the sign for reporting

    % ---- exposure: fraction of the segment above zcut -------------------
    expose = local_exposedFraction(P1(3), P2(3), zcut);

    n1  = N1(ii,:);  n2 = N2(ii,:);
    bis = n1 + n2;
    bis = bis / norm(bis);

    amp = cube.edgeAmp(ii) * cube.edgeLenW(ii) * expose^cube.exposeExp;

    E(ii).family   = cube.family(ii);
    E(ii).axisId   = cube.edgeAxis(ii);
    E(ii).axis     = t;
    E(ii).n1       = n1;
    E(ii).n2       = n2;
    E(ii).bisector = bis;
    E(ii).expose   = expose;
    E(ii).amp      = amp;

    % ---- diagnostics ----------------------------------------------------
    E(ii).tilt = asind(min(abs(t(3)),1));
    % The stripe is the great circle perpendicular to t; where it crosses the
    % DRP it runs perpendicular to t's horizontal projection.
    E(ii).stripeAz = mod(atan2d(t(2), t(1)) + 90, 180);
    % A surface whose normal is n reflects into the DRP cell
    %   theta = 2*(elevation(n) - 45),  azimuth = atan2(n_y, n_x)
    % A bevel facing below the horizon cannot reflect into the camera at all,
    % so it has no DRP theta - report NaN rather than a meaningless negative.
    E(ii).arcCentreAz = mod(atan2d(bis(2), bis(1)), 360);
    if bis(3) > 0
        E(ii).arcCentreTh = 2 * (asind(min(bis(3),1)) - 45);
    else
        E(ii).arcCentreTh = NaN;
    end
    E(ii).zMid        = 0.5 * (P1(3) + P2(3));

    E(ii).cAxis = cS;
    E(ii).cAz   = cAz;
    E(ii).Phi   = Phi;
end
end

% =========================================================================
function f = local_exposedFraction(z1, z2, zcut)
% Fraction of the segment z1 -> z2 that lies strictly above zcut.  z is linear
% in the segment parameter, so this is exact.
if z1 == z2
    f = double(z1 > zcut);
    return
end
s0 = (zcut - z1) / (z2 - z1);          % parameter where the edge crosses zcut
if z2 > z1
    f = 1 - s0;                        % exposed for s > s0
else
    f = s0;                            % exposed for s < s0
end
f = min(max(f,0),1);
end

% =========================================================================
function e = local_blank()
e = struct( ...
    'family',       "", ...
    'axisId',       0, ...
    'axis',         [NaN NaN NaN], ...
    'n1',           [NaN NaN NaN], ...
    'n2',           [NaN NaN NaN], ...
    'bisector',     [NaN NaN NaN], ...
    'expose',       0, ...
    'amp',          0, ...
    'tilt',         NaN, ...
    'stripeAz',     NaN, ...
    'arcCentreAz',  NaN, ...
    'arcCentreTh',  NaN, ...
    'zMid',         NaN, ...
    'cAxis',        [NaN NaN NaN], ...
    'cAz',          NaN, ...
    'Phi',          NaN);
end
