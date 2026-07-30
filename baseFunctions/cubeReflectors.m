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
% THE SURFACE-TRACE APPROXIMATION ----------------------------------------
% What reflects is not the edge of the buried box but the CORRUGATION CREST
% that edge cuts in the polished surface, and a surface trace is horizontal by
% construction.  So the crest direction used for the optics is the edge axis
% projected into the surface plane, not the tilted edge itself.
%
% This matters because the stripe is the great circle perpendicular to the
% crest: a HORIZONTAL crest gives a great circle containing +Z, i.e. a stripe
% that is a DIAMETER THROUGH THE DRP CENTRE - the defining observation of
% section 2.3.2.  A tilted crest would give an off-centre chord, which is not
% what a DRP of an etched surface shows.
%
% The bevel arc is projected with it: .arcN1 and .arcN2 are the two face
% normals brought into the plane perpendicular to the crest, so the lit arc
% still runs from one face round to the other, and its half-angle .arcWedge is
% no longer exactly 90 deg once the edge was tilted.
%
%   E  12x1 struct array, one per edge:
%     .family       "broadEnd" (axis || c) | "basalEnd" | "basalBroad"
%     .axisId       1 = || e_c, 2 = || e_b, 3 = || e_p
%     .axis         1x3 unit TRUE edge direction, SAMPLE frame, z >= 0
%     .crest        1x3 unit HORIZONTAL surface trace of that edge - this is
%                   what DRPsim_cube uses, and why every stripe is a diameter
%     .n1 .n2       1x3 outward normals of the two faces the bevel joins
%     .arcN1 .arcN2 1x3 those normals projected into the plane perpendicular
%                   to .crest: the two ends of the lit arc
%     .arcWedge     deg, angle between them (90 for an untilted edge)
%     .bisector     1x3 normalise(arcN1+arcN2): the direction the bevel faces
%     .expose       0..1 fraction of the edge above the etch cut
%     .amp          edgeAmp * lengthWeight * expose^exposeExp
%     .tilt         deg, how far the TRUE edge tilts out of the surface plane.
%                   Now a pure diagnostic: it says how much projecting onto
%                   the surface changed the crest.  At tilt = 90 the edge is
%                   vertical, its trace is a point, and the edge is silenced.
%     .stripeAz     deg mod 180, azimuth of the stripe (= crest azimuth + 90)
%     .arcCentreAz  deg 0..360, azimuth where the bevel faces
%     .arcCentreTh  deg, DRP theta of the same, = 2*(elevation(bisector) - 45),
%                   or NaN if the bevel faces below the horizon.  This is where
%                   the bevel puts its MAXIMUM along the stripe, i.e. the
%                   asymmetry a real rounded edge introduces.
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

    % ---- surface trace: the crest is the edge projected into the surface -
    % A crest cut into a polished surface is horizontal, so this is what the
    % optics must use.  It is also what forces the stripe through the DRP
    % centre.  A vertical edge has no trace and cannot make a stripe.
    cr = [t(1) t(2) 0];
    crN = norm(cr);
    if crN < 1e-9
        continue                        % edge is vertical: no surface trace
    end
    cr = cr / crN;

    % the bevel arc travels with it: bring both face normals into the plane
    % perpendicular to the crest
    a1 = local_projectPerp(n1, cr);
    a2 = local_projectPerp(n2, cr);
    if isempty(a1) || isempty(a2)
        continue                        % a face normal is along the crest
    end
    wedge = acosd(min(max(dot(a1,a2),-1),1));
    bis   = a1 + a2;
    if norm(bis) < 1e-9
        continue                        % the two faces face opposite ways
    end
    bis = bis / norm(bis);

    amp = cube.edgeAmp(ii) * cube.edgeLenW(ii) * expose^cube.exposeExp;

    E(ii).family   = cube.family(ii);
    E(ii).axisId   = cube.edgeAxis(ii);
    E(ii).axis     = t;
    E(ii).crest    = cr;
    E(ii).n1       = n1;
    E(ii).n2       = n2;
    E(ii).arcN1    = a1;
    E(ii).arcN2    = a2;
    E(ii).arcWedge = wedge;
    E(ii).bisector = bis;
    E(ii).expose   = expose;
    E(ii).amp      = amp;

    % ---- diagnostics ----------------------------------------------------
    E(ii).tilt = asind(min(abs(t(3)),1));
    % The stripe is the great circle perpendicular to the crest.  The crest is
    % horizontal, so that circle contains +Z and the stripe is a diameter.
    E(ii).stripeAz = mod(atan2d(cr(2), cr(1)) + 90, 180);
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
function u = local_projectPerp(v, ax)
% Component of v perpendicular to the unit vector ax, renormalised.  Empty if
% v is parallel to ax (nothing left to normalise).
u = v - dot(v,ax)*ax;
n = norm(u);
if n < 1e-9
    u = [];
else
    u = u / n;
end
end

% =========================================================================
function e = local_blank()
e = struct( ...
    'family',       "", ...
    'axisId',       0, ...
    'axis',         [NaN NaN NaN], ...
    'crest',        [NaN NaN NaN], ...
    'n1',           [NaN NaN NaN], ...
    'n2',           [NaN NaN NaN], ...
    'arcN1',        [NaN NaN NaN], ...
    'arcN2',        [NaN NaN NaN], ...
    'arcWedge',     NaN, ...
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
