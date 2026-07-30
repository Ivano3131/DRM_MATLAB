function [simDRP, E] = DRPsim_cube(eu1, eu2, eu3, exp_para)
% DRPsim_cube  Simulate the DRP of an etched Ti64 alpha lath modelled as a box.
%
%   simDRP      = DRPsim_cube(eu1,eu2,eu3,exp_para)
%   [simDRP, E] = DRPsim_cube(...)
%
% Paints the exposed beveled edges returned by cubeReflectors onto the same
% th_num x ph_num grid, with the same probe-direction convention, as DRPsim /
% DRPsim_hcp / DRPsim_Ti64 - so DRPdisp, drp_loader, igray2drp and the
% dictionary code all work unchanged.
%
% HOW AN EDGE IS PAINTED --------------------------------------------------
% Every DRP cell (theta,phi) stands for the microfacet normal that would send
% the illumination specularly into the fixed overhead camera:
%
%       v(theta,phi) = [cos(phi)cos(45+theta/2), sin(phi)cos(45+theta/2), sin(45+theta/2)]
%
% A rounded edge is a piece of cylinder about its CREST t, so its surface
% normals all lie on the great circle perpendicular to t.  The crest is the
% edge's HORIZONTAL surface trace (cubeReflectors projects it), because what
% reflects is the corrugation cut into the polished surface, not the tilted
% edge of the buried box.  Being horizontal, that great circle contains +Z, so
% EVERY STRIPE IS A DIAMETER THROUGH THE DRP CENTRE - the defining observation
% of section 2.3.2.
%
%   ACROSS the stripe - the angular distance of a probe from that great
%   circle is asind(v.t), giving a Cauchy ridge of half-width stripeWidth:
%
%       across = amp / (1 + (asind(v.t)/stripeWidth)^2)
%
%   ALONG the stripe - the bevel is only an ARC of that cylinder: rounding the
%   edge sweeps the normal from one face round to the other and no further.
%   Projecting the probe onto the plane perpendicular to the crest,
%
%       w   = normalise(v - (v.t)t)
%       psi = angle of w from f.arcN1 towards f.arcN2
%
%   the edge lights psi within [0, arcWedge] widened by bevelSpan, with a
%   raised-cosine shoulder of width bevelSoft.  arcWedge is 90 deg for an
%   untilted edge and less once the edge has been projected onto the surface.
%   This is why the stripe is BRIGHTER ON ONE SIDE of the DRP centre: the lit
%   arc is centred on the direction the bevel faces, not on +Z.
%
%   Finally a (v.Z)^stripeTaper microfacet foreshortening factor, the same
%   term as in DRPsim_Ti64.
%
% Edges combine with max() (brightest reflector wins) unless
% exp_para.cube.combine == "sum".  The result is min-subtracted and
% peak-normalised to [0,255], so it is the RATIO of the two stripe families
% that carries the orientation information, not their absolute values.
%
% INPUTS ------------------------------------------------------------------
%   eu1,eu2,eu3  Bunge Euler angles [phi1 Phi phi2], degrees.
%   exp_para     needs th_min, th_max, th_num, ph_num and .cube
%                (cubeGeometryTi64).
%
% OUTPUTS -----------------------------------------------------------------
%   simDRP  th_num x ph_num uint8; rows = theta ascending from th_min,
%           columns = azimuth 0 : 360/ph_num : 360-360/ph_num.
%   E       the edge list used, for verification (see cubeReflectors).
%
% See also cubeGeometryTi64, cubeReflectors, checkCubeMechanism,
%          makeDRPdic_cube, matchDRPcube.
% -------------------------------------------------------------------------
arguments
    eu1 double
    eu2 double
    eu3 double
    exp_para struct
end

if ~isfield(exp_para,'cube') || isempty(exp_para.cube)
    error('DRPsim_cube:notSetup', ...
        'exp_para.cube is missing. Run exp_para = cubeGeometryTi64(exp_para) first.');
end
cube = exp_para.cube;

th_min = exp_para.th_min;
th_max = exp_para.th_max;
th_num = exp_para.th_num;
ph_num = exp_para.ph_num;

% ---- DRP sampling grid --------------------------------------------------
% The azimuth step MUST be 360/ph_num: matchDRPcube recovers phi1 by shifting
% whole columns, which is only a rotation if the columns close the circle.
th_step = (th_max - th_min) / (th_num - 1);
ph_step = 360 / ph_num;
th_DRP  = repmat((th_min : th_step : th_max).', 1, ph_num);
ph_DRP  = repmat(0 : ph_step : 360-ph_step, th_num, 1);

% ---- probe direction of every cell: the specular microfacet normal ------
elev = 45 + th_DRP/2;
vx   = cosd(ph_DRP) .* cosd(elev);
vy   = sind(ph_DRP) .* cosd(elev);
vz   = sind(elev);

if cube.stripeTaper > 0
    taper = max(vz,0) .^ cube.stripeTaper;
else
    taper = 1;
end

% ---- the physics: one edge list for this orientation --------------------
E = cubeReflectors(eu1, eu2, eu3, cube);

simDRP = zeros(th_num, ph_num);
useSum = cube.combine == "sum";

for ii = 1:numel(E)
    f = E(ii);
    if ~(f.amp > 0), continue; end

    % The CREST is the edge's horizontal surface trace, not the tilted edge
    % itself (see cubeReflectors).  Being horizontal, its great circle
    % contains +Z, so every stripe is a diameter through the DRP centre.
    t = f.crest;

    % (a) across the stripe: distance from the great circle with pole t
    vt   = t(1)*vx + t(2)*vy + t(3)*vz;
    dist = asind(min(max(vt,-1),1));
    temp = f.amp ./ (1 + (dist ./ cube.stripeWidth).^2);

    % (b) along the stripe: the bevel gate
    temp = temp .* local_bevelGate(f, t, vt, vx, vy, vz, ...
                                   cube.bevelSpan, cube.bevelSoft);

    % (c) microfacet foreshortening
    temp = temp .* taper;

    if useSum
        simDRP = simDRP + temp;
    else
        simDRP = max(simDRP, temp);
    end
end

% ---- normalise to [0,255] ----------------------------------------------
simDRP = simDRP - min(simDRP(:));
peak   = max(simDRP(:));
if peak > 0
    simDRP = simDRP / peak;
end
simDRP = uint8(simDRP * 255);
end

% =========================================================================
function g = local_bevelGate(f, t, vt, vx, vy, vz, span, soft)
% Fraction of the edge's rounding that reflects into each DRP cell.
%
% The bevel's normals run from arcN1 to arcN2 along the great circle
% perpendicular to the crest, sweeping f.arcWedge degrees (90 for an untilted
% edge, less once the edge was projected onto the surface).  Project the probe
% into that plane and read off its position psi; the edge lights psi within
% (arcWedge + span)/2 of the wedge centre, fading to zero over a further
% `soft` degrees.  A span large enough to reach 180 opens the gate to the
% whole circle, restoring the idealised symmetric stripe of a semi-cylinder.
wedge    = f.arcWedge;
halfSpan = wedge/2 + span/2;
if isnan(wedge) || halfSpan >= 180
    g = 1;
    return
end

% probe projected onto the plane perpendicular to the crest, read in an
% orthonormal basis of that plane built on arcN1
wx = vx - vt*t(1);
wy = vy - vt*t(2);
wz = vz - vt*t(3);

e1 = f.arcN1;
e2 = f.arcN2 - dot(f.arcN2,e1)*e1;        % arcN1 and arcN2 need not be perpendicular
n2n = norm(e2);
if n2n < 1e-9
    g = 1;
    return
end
e2 = e2 / n2n;

c1 = wx*e1(1) + wy*e1(2) + wz*e1(3);
c2 = wx*e2(1) + wy*e2(2) + wz*e2(3);

psi  = atan2d(c2, c1);                              % 0 at arcN1, arcWedge at arcN2
dpsi = abs(mod(psi - wedge/2 + 180, 360) - 180);    % wrapped distance from the wedge centre

if soft > 0
    x = min(max((dpsi - halfSpan) ./ soft, 0), 1);
    g = 0.5 * (1 + cos(pi*x));
else
    g = double(dpsi <= halfSpan);
end

% where the probe is parallel to t the projection is undefined, but `across`
% is already ~0 there, so leaving the gate open costs nothing
g(c1 == 0 & c2 == 0) = 1;
end
