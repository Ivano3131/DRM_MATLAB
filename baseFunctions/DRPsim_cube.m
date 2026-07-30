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
% A rounded edge with axis t is a piece of cylinder, so its surface normals
% all lie on the great circle perpendicular to t.
%
%   ACROSS the stripe - the angular distance of a probe from that great
%   circle is asind(v.t), giving a Cauchy ridge of half-width stripeWidth:
%
%       across = amp / (1 + (asind(v.t)/stripeWidth)^2)
%
%   ALONG the stripe - the bevel is only a QUARTER of that cylinder: rounding
%   the edge sweeps the normal from face n1 round to face n2 and no further.
%   Projecting the probe onto the plane perpendicular to t,
%
%       w   = normalise(v - (v.t)t)
%       psi = atan2d(w.n2, w.n1)          0 at face n1, 90 at face n2
%
%   and the edge only lights psi within [0,90] widened by bevelSpan, with a
%   raised-cosine shoulder of width bevelSoft.  This is why the stripe is
%   BRIGHTER ON ONE SIDE of the DRP centre: the lit arc is centred on the
%   direction the bevel faces, not on +Z.
%
%   Note the parameterisation is in the plane perpendicular to t, NOT in the
%   {Z, u} plane DRPsim_Ti64 uses.  That matters here because a box edge is
%   free to tilt out of the surface, and a tilted crest's great circle does
%   not pass through the DRP centre.
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

    t = f.axis;

    % (a) across the stripe: distance from the great circle with pole t
    vt   = t(1)*vx + t(2)*vy + t(3)*vz;
    dist = asind(min(max(vt,-1),1));
    temp = f.amp ./ (1 + (dist ./ cube.stripeWidth).^2);

    % (b) along the stripe: the quarter-cylinder bevel gate
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
% The bevel's normals run from n1 to n2 along the great circle perpendicular
% to t, a 90 deg quarter.  Project the probe into that plane and read off its
% position psi; the edge lights psi within (90 + span)/2 of the wedge centre
% (psi = 45), fading to zero over a further `soft` degrees.  A span >= 270
% opens the gate to the whole circle, restoring the idealised symmetric
% stripe of a full semi-cylinder.
halfSpan = 45 + span/2;
if halfSpan >= 180
    g = 1;
    return
end

% probe projected onto the plane perpendicular to t, then read in the {n1,n2}
% basis (n1 and n2 are orthonormal there, being two box axes)
wx = vx - vt*t(1);
wy = vy - vt*t(2);
wz = vz - vt*t(3);

n1 = f.n1;  n2 = f.n2;
c1 = wx*n1(1) + wy*n1(2) + wz*n1(3);
c2 = wx*n2(1) + wy*n2(2) + wz*n2(3);

psi  = atan2d(c2, c1);                              % 0 at n1, 90 at n2
dpsi = abs(mod(psi - 45 + 180, 360) - 180);         % wrapped distance from the wedge centre

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
