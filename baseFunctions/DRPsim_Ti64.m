function [simDRP, R] = DRPsim_Ti64(eu1,eu2,eu3,exp_para)
% DRPsim_Ti64  Simulate the DRP of an etched wrought-Ti64 alpha lath.
%
%   simDRP      = DRPsim_Ti64(eu1,eu2,eu3,exp_para)
%   [simDRP, R] = DRPsim_Ti64(...)
%
% Renders the reflecting features returned by lathReflectors into a
% th_num x ph_num uint8 DRP on the same grid, with the same probe-direction
% convention, as DRPsim / DRPsim_hcp, so the existing dictionary, indexing
% engine and display code work unchanged.
%
% HOW EACH FEATURE IS PAINTED --------------------------------------------
% Every DRP cell (theta,phi) stands for the microfacet normal that would send
% the illumination specularly into the fixed overhead camera:
%
%       v(theta,phi) = thph2vec(45 + theta/2, phi)
%
%   "ridge" features (rounded crests)
%       A cylinder with axis t reflects the camera along the whole great
%       circle perpendicular to t.  The angular distance of a probe from that
%       great circle is asind(v . t), so ACROSS the stripe the profile is
%           amp / (1 + (asind(v.t)/stripeWidth)^2)
%       a Cauchy ridge of half-width stripeWidth.  Because t is horizontal,
%       the great circle contains +Z and the stripe is a diameter through the
%       DRP centre - the defining observation of section 2.3.2.
%
%       ALONG the stripe the profile is NOT uniform.  An etched crest is a
%       BEVEL joining two facets, not a full semi-cylinder, so its normals
%       sweep only the wedge between those facets and it lights only an ARC
%       of the great circle.  The arc is centred on the direction the crest
%       faces (lathReflectors computes it per feature: the prism wedge
%       bisector for the terrace ridge, the basal terrace normal for the
%       basal edge), which is tilted away from +Z by the facet angles - so
%       the stripe comes out BRIGHTER ON ONE SIDE of the DRP centre.  Two
%       further factors shape the profile: the raised-cosine shoulders of
%       that arc (bevelSoft) and the (n.Z)^stripeTaper microfacet
%       foreshortening.  Set bevelSpan >= 360 to recover the idealised
%       symmetric stripe.
%
%   "peak" features (flat plateaus / pyramidal apices)
%       A plane mirror with normal n is bright only where the probe equals n,
%       giving a point lobe
%           amp / (1 + (acosd(v.n)/peakWidth)^2)
%       which lands at theta = 2*(elevation(n) - 45).
%
% Features are combined with max() by default (brightest reflector wins, the
% same rule as DRPsim/DRPsim_hcp) or by incoherent addition when
% exp_para.lath.combine == "sum".  The finished pattern is rescaled to
% [0,255], so it is the RATIO of the two stripe amplitudes that carries the
% orientation information, not their absolute values.
%
% INPUTS ------------------------------------------------------------------
%   eu1,eu2,eu3  Bunge Euler angles [phi1 Phi phi2] in degrees.
%   exp_para     needs the DRP grid (th_min, th_max, th_num, ph_num) and
%                exp_para.lath from lathGeometryTi64.
%
% OUTPUTS -----------------------------------------------------------------
%   simDRP  th_num x ph_num uint8, rows = theta ascending from th_min,
%           columns = azimuth 0 : 360/ph_num : 360-360/ph_num.
%   R       the reflector list used, for verification (see lathReflectors).
%
% See also lathGeometryTi64, lathReflectors, checkTi64Mechanism,
%          makeDRPdic_Ti64, check_indexing_result_Ti64.
% -------------------------------------------------------------------------
arguments
    eu1 double
    eu2 double
    eu3 double
    exp_para struct
end

if ~isfield(exp_para,'lath') || isempty(exp_para.lath)
    error('DRPsim_Ti64:notSetup', ...
        'exp_para.lath is missing. Run exp_para = lathGeometryTi64(exp_para) first.');
end
lath = exp_para.lath;

th_min = exp_para.th_min;
th_max = exp_para.th_max;
th_num = exp_para.th_num;
ph_num = exp_para.ph_num;

% ---- DRP sampling grid --------------------------------------------------
% The azimuth step MUST be 360/ph_num: the dictionary and IndexingEngine
% align patterns by circshifting whole columns.
th_step = (th_max - th_min) / (th_num - 1);
ph_step = 360 / ph_num;
th_DRP  = repmat((th_min : th_step : th_max).', 1, ph_num);
ph_DRP  = repmat(0 : ph_step : 360-ph_step, th_num, 1);

% ---- probe direction of every cell: the specular microfacet normal ------
elev = 45 + th_DRP/2;
vx   = cosd(ph_DRP) .* cosd(elev);
vy   = sind(ph_DRP) .* cosd(elev);
vz   = sind(elev);

cauchy = @(amp,width,dist) amp ./ (1 + (dist./width).^2);

% ---- optional taper: fade stripes toward the DRP edge (low theta) -------
if lath.stripeTaper > 0
    taper = max(vz,0) .^ lath.stripeTaper;
else
    taper = 1;
end

% ---- the physics: one reflector list for this orientation ---------------
R = lathReflectors(eu1, eu2, eu3, lath);

simDRP  = zeros(th_num, ph_num);
useSum  = lath.combine == "sum";

for ii = 1:numel(R)
    f = R(ii);
    if ~(f.amp > 0), continue; end

    switch f.kind
        case "ridge"
            % (a) ACROSS the stripe: distance of each probe from the great
            %     circle whose pole is the (horizontal) crest direction.
            t    = f.crest;
            dot3 = t(1)*vx + t(2)*vy + t(3)*vz;
            dist = asind(max(min(dot3,1),-1));
            temp = cauchy(f.amp, lath.stripeWidth, dist) .* taper;

            % (b) ALONG the stripe: the bevel gate.  A real crest joins two
            %     facets, so its normals sweep only the wedge between them,
            %     not the full semi-cylinder; it therefore lights only an arc
            %     of the great circle, centred on the direction it faces
            %     (f.arcCentre).  Without this the stripe would be uniform
            %     along the circle and hence perfectly symmetric about the
            %     DRP centre, which a beveled crest is not.
            temp = temp .* local_bevelGate(f, vx, vy, vz, lath.bevelSoft);

        case "peak"
            n    = f.normal;
            dot3 = n(1)*vx + n(2)*vy + n(3)*vz;
            dist = acosd(max(min(dot3,1),-1));
            temp = cauchy(f.amp, lath.peakWidth, dist);

        otherwise
            continue
    end

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
function g = local_bevelGate(f, vx, vy, vz, soft)
% Fraction of the crest's rounding that actually reflects into each DRP cell.
%
% Position along the stripe is the angle psi in the {Z, arcU} plane:
%       psi = atan2(v . arcU, v . Z)
% and the rounded region covers psi within f.arcSpan/2 of f.arcCentre, fading
% to zero over a further `soft` degrees (raised cosine, so nothing is
% hard-cut - real bevels have rounded shoulders).  Both the centre and the
% span come from the facets the crest joins (see lathReflectors); a span
% >= 360 restores the full semi-cylinder, i.e. a stripe that is uniform and
% symmetric about the DRP centre.
span = f.arcSpan;
if isnan(span) || span <= 0 || span >= 360 || isnan(f.arcCentre) || any(isnan(f.arcU))
    g = 1;
    return
end
u    = f.arcU;
psi  = atan2d(u(1)*vx + u(2)*vy + u(3)*vz, vz);
dpsi = abs(mod(psi - f.arcCentre + 180, 360) - 180);   % wrapped separation
if soft > 0
    x = min(max((dpsi - span/2) ./ soft, 0), 1);
    g = 0.5 * (1 + cos(pi*x));
else
    g = double(dpsi <= span/2);
end
end
