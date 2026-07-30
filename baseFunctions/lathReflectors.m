function R = lathReflectors(eu1,eu2,eu3,lath)
% lathReflectors  Specimen-frame reflecting features of an etched Ti64 lath.
%
%   R = lathReflectors(eu1,eu2,eu3,lath)
%
% This is where the physics of section 2.3.2 lives.  Given one alpha
% orientation (Bunge [phi1 Phi phi2] in degrees) and the model built by
% lathGeometryTi64, it returns every reflecting feature with ALL of its
% intermediate quantities exposed, so the mechanism can be verified number by
% number without ever rendering a DRP.  DRPsim_Ti64 simply paints the lobes
% that this function describes.
%
% NO MTEX is used here (the crystal-frame vectors were precomputed), so the
% function is safe inside parfor.
%
% GEOMETRY RECAP ----------------------------------------------------------
%   c_S            the c-axis in the specimen frame, forced to point up.
%   Phi            angle between c_S and the sample normal Z.  For the Bunge
%                  convention this is the second Euler angle (folded to
%                  [0,90] because c is an axis, not a vector).
%   c_h            in-surface projection of c_S; its azimuth cAz is the
%                  "c-axis direction" seen in the DRP.  cAz = phi1 - 90 deg.
%
%   A rounded (semi-cylindrical) crest along a horizontal unit axis t
%   reflects the fixed overhead camera along the whole great circle
%   perpendicular to t.  In the DRP that great circle is the DIAMETER at
%   azimuth (azimuth of t) + 90 deg.  Hence:
%       crest || c        ->  stripe PERPENDICULAR to the c-axis direction
%       crest perp. to c  ->  stripe PARALLEL      to the c-axis direction
%
%   A crest axis that is NOT horizontal cannot produce a corrugation whose
%   crest line lies in the surface; the etched feature degenerates from a
%   semi-cylinder into a semi-spherical (pyramidal) peak.  This is captured
%   by two separate numbers per feature:
%       tilt   angle of the crystallographic crest axis out of the surface
%       coinc  = cos(tilt)^coincExp, in [0,1]; 1 = perfect semi-cylinder
%   The observable crest direction is always the in-surface PROJECTION of the
%   axis, so every stripe still passes through the DRP centre, as observed.
%
% FIELDS OF R (struct array, one element per feature) ---------------------
%   name      "terraceRidge" | "basalEdge" | "basalPlateau" | "prismPlateau"
%   kind      "ridge" (great-circle stripe) | "peak" (point maximum)
%   variant   index into lath.variantRot (1 = the indexed alpha orientation)
%   axis      1 x 3  crystallographic crest axis in the specimen frame (ridge)
%   crest     1 x 3  unit in-surface crest direction  = normalised proj. of axis
%   crestAz   azimuth of crest, deg
%   stripeAz  azimuth of the DRP stripe = crestAz + 90, folded to [0,180)
%   tilt      angle of axis out of the sample surface, deg
%   coinc     cos(tilt)^coincExp, the semi-cylinder / pyramid split
%   arcU      1 x 3  in-surface basis vector of the stripe's great circle;
%             positions along the stripe are angles in the {Z, arcU} plane
%   arcCentre where along the stripe the BEVEL faces, deg from the DRP centre
%             (+ve toward stripeAz, -ve toward stripeAz + 180).  This is what
%             makes the stripe asymmetric: a real crest is a bevel joining two
%             facets, so it only lights the arc it faces, not the whole circle
%   arcCentreTh / arcCentreAz  the same position as DRP (theta, azimuth), i.e.
%             where the stripe should be brightest.  Compare against the data
%   arcWedge  the crystallographic wedge the surface turns through at this
%             crest, deg (angle between the two facets it joins).  0 when the
%             second facet is not identified
%   arcSpan   the arc actually lit = arcWedge + bevelSpan, the second term
%             being the extra rounding the etch adds.  Half of it either side
%             of arcCentre at full intensity, then a bevelSoft shoulder
%   edgeKind  basalEdge only: "prism" if this rim is bounded by a prismatic
%             side face, "lil" if by a side face along the lattice invariant
%             line.  A rectangular lath has both
%   rim       basalEdge only: +1 / -1, the two parallel rims of the basal
%             rectangle (their bevels face opposite ways, so together they
%             give the parallel-to-c stripe its two lobes)
%   normal    1 x 3  specimen-frame plane normal (peak features)
%   peakTh    DRP theta at which the specular point appears, deg (peaks);
%             = 2*(elevation of normal - 45).  NaN when unreachable.
%   peakAz    DRP azimuth of the specular point, deg (peaks)
%   expose    projected-area factor of the parent facet, in [0,1]
%   shadow    low-elevation shadowing factor, in [0,1] (peaks)
%   amp       final relative amplitude fed to DRPsim_Ti64
%   note      short human-readable statement of what this feature is
%
% See also lathGeometryTi64, DRPsim_Ti64, checkTi64Mechanism.
% -------------------------------------------------------------------------
arguments
    eu1 double
    eu2 double
    eu3 double
    lath struct
end

Z   = [0 0 1];
tol = 1e-9;

R = repmat(local_blank(), 0, 1);  % empty 0 x 1 struct array with the right fields
k = 0;

for v = 1:numel(lath.variantRot)

    % ---------------------------------------------------------------------
    % Crystal-frame vectors of this alpha variant, expressed in the frame of
    % the indexed orientation.  Rv = eye(3) for v = 1 (the grain itself);
    % for v > 1 it is the Burgers-OR misorientation to a sibling variant.
    % Row-vector convention: u_row * Rv.'  ==  (Rv * u_col).'
    % ---------------------------------------------------------------------
    Rv = lath.variantRot{v};
    wv = lath.variantW(v);

    cX     = lath.cAxis        * Rv.';
    ledX   = lath.ledgeDirs    * Rv.';
    prismX = lath.prismNormals * Rv.';

    % ---- crystal -> specimen (active Bunge, same as the rest of the repo)
    cS = EulerRotate(cX, eu1, eu2, eu3);
    cS = cS / norm(cS);
    if cS(3) < 0, cS = -cS; end               % c is an axis: keep the upper one

    Phi = acosd(max(min(cS(3),1),-1));        % angle(c, sample normal), [0,90]

    cH  = cS - cS(3)*Z;                       % in-surface projection of c
    nH  = norm(cH);
    if nH > tol
        cHu = cH / nH;
        cAz = atan2d(cHu(2), cHu(1));         % = phi1 - 90 deg
    else
        cHu = [NaN NaN NaN];                  % c || Z: no in-surface c direction
        cAz = NaN;
    end

    % ---- how much basal / prism area the sample surface exposes ----------
    % A plane with normal n presents a projected area proportional to |n.Z|.
    % Basal normal is c            -> cos(Phi).
    % Prism normals are perp. to c -> at most sin(Phi).
    exposeBasal = cosd(Phi) ^ lath.exposeExp;
    exposePrism = sind(Phi) ^ lath.exposeExp;

    % ---- prism plane normals in the specimen frame, upward-facing --------
    % Needed by BOTH the terrace ridge (they bound its bevel) and the prism
    % plateau peaks, so they are built once here.
    prismS = normr(EulerRotate(prismX, eu1, eu2, eu3));
    flip   = prismS(:,3) < 0;
    prismS(flip,:) = -prismS(flip,:);

    % =====================================================================
    % FEATURE 1  "terraceRidge":  crest || c  ->  stripe PERPENDICULAR to c
    % ---------------------------------------------------------------------
    % The ridge that joins the two exposed prism facets of the lath (Fig
    % 2.15, point 1: it connects (10-10) and (11-20)).  All prism planes
    % contain c, so every prism-prism intersection line is parallel to c.
    %   axis  = c_S
    %   tilt  = 90 - Phi  (c lies in the surface at Phi = 90 -> best ridge)
    %   amp  ~ sin(Phi)^exposeExp   (prism facets must be exposed)
    %        * cos(tilt)^coincExp   (ridge line must lie in the surface)
    % Both factors reduce to sin(Phi) here: this feature switches on as the
    % c-axis rotates into the surface, which is why it dominates at point 1
    % (Phi = 87.9 deg) and is only a faint second stripe at point 4
    % (Phi = 37.8 deg).  Set coincExp = 0 to keep the exposure term alone.
    %
    % BEVEL: the crest is not a full semi-cylinder.  Its normals sweep only
    % the wedge between the two prism facets it joins, so it lights only an
    % ARC of the great circle.  The arc is centred on the wedge bisector of
    % the two most exposed prism facets - i.e. on the direction the crest
    % actually faces - which is generally tilted away from +Z, so the stripe
    % comes out brighter on one side of the DRP centre than the other.
    % =====================================================================
    tiltT  = asind(min(abs(cS(3)),1));                % == 90 - Phi
    coincT = cosd(tiltT) ^ lath.coincExp;
    ampT   = lath.ampTerraceRidge * wv * exposePrism * coincT;
    if nH > tol && ampT > 0
        uT = local_arcBasis(cHu);
        % The two most upward-facing prism facets are the pair the crest
        % joins.  Their bisector is the direction the crest FACES (arc
        % centre) and the angle between them is the wedge the surface must
        % turn through (arc width before etch rounding).
        [~,ordT] = sort(prismS(:,3),'descend');
        n1 = prismS(ordT(1),:);
        if size(prismS,1) >= 2
            n2 = prismS(ordT(2),:);
        else
            n2 = n1;
        end
        bis   = n1 + n2;
        if norm(bis) > tol, bis = bis / norm(bis); else, bis = Z; end
        wedgeT = acosd(max(min(dot(n1,n2),1),-1));
        % The two prism facets generally STRADDLE the surface normal, so their
        % bisector sits close to +Z and this crest faces almost straight up ->
        % a nearly symmetric stripe.  The real crest is not symmetric: one of
        % its flanks is the etch-resistant beta lamella (Fig 2.15 shows a
        % clearly asymmetric sawtooth, and Fig 2.14b puts the lath broad face
        % ~14 deg off the ledge terrace normal).  That flank cannot be derived
        % from the alpha orientation alone without committing to a parent beta
        % variant, so it is left as a signed fitted tilt.
        psiT = local_psi(bis, uT) + lath.bevelBias;

        f = local_blank();
        f.name     = "terraceRidge";
        f.kind     = "ridge";
        f.variant  = v;
        f.axis     = cS;
        f.crest    = cHu;
        f.crestAz  = cAz;
        f.stripeAz = mod(cAz + 90, 180);
        f.tilt     = tiltT;
        f.coinc    = coincT;
        f.expose   = exposePrism;
        f.amp      = ampT;
        f.arcU     = uT;
        f.arcCentre   = psiT;
        f.arcCentreTh = local_arcTheta(f.arcCentre);
        f.arcCentreAz = local_arcAz(f.arcCentre, f.arcU);
        f.arcWedge    = wedgeT;
        f.arcSpan     = wedgeT + lath.bevelSpan;
        f.note     = "prism-prism ridge, axis || c  -> stripe PERP to c";
        k = k + 1;  R(k,1) = f;
    end

    % =====================================================================
    % FEATURE 2  "basalEdge":  crest perp. to c  ->  stripe PARALLEL to c
    % ---------------------------------------------------------------------
    % The rounded edge between the chemically resistant beta matrix and the
    % slow-etching basal terrace (Fig 2.15, point 4).  It lies IN the basal
    % plane, so it is perpendicular to c, and its stripe is parallel to the
    % c-axis direction.
    %   amp ~ cos(Phi)^exposeExp   (basal must be exposed -> small Phi)
    %       * cos(tilt)^coincExp   (edge line must lie in the surface)
    %
    % WHICH basal direction?  An alpha lath is a RECTANGULAR platelet with 6
    % faces (2 basal broad faces, 2 prismatic side faces, 2 side faces along
    % the lattice invariant line - Fig 2.14), so its basal face is an
    % approximate rectangle with exactly TWO rim directions: one prismatic,
    % one along the lattice invariant line.  A given grain still has several
    % candidate laths, because the parent beta is not determined by the alpha
    % orientation alone, so ledgeSelect chooses among them:
    %   "best"  keep one lath - in the rectangular habit BOTH of its rim
    %           directions are emitted, in the hexagonal habit just the single
    %           best-coinciding edge.  This is the phi2-sensitive branch:
    %           rotating the cell about c moves the chosen rim off the exact
    %           basal trace, which dims the stripe and swings its azimuth off
    %           exactly-parallel-to-c - one mechanism for the reported
    %           non-orthogonal stripes.
    %   "trace" use the exact basal surface trace c x Z: always horizontal,
    %           stripe exactly parallel to the c-axis, no phi2 dependence.
    %   "all"   every candidate (diagnostic).
    % =====================================================================
    if lath.ledgeSelect == "trace"
        tBasal = cross(cS, Z);                        % basal trace, horizontal
        if norm(tBasal) > tol
            edgeAxes = tBasal / norm(tBasal);
            edgeW    = 1;  edgeKind = 1;  edgePair = 1;
        else
            edgeAxes = zeros(0,3);
            edgeW = zeros(0,1);  edgeKind = zeros(0,1);  edgePair = zeros(0,1);
        end
    else
        edgeAxes = normr(EulerRotate(ledX, eu1, eu2, eu3));
        edgeW    = lath.ledgeW;
        edgeKind = lath.ledgeKind;
        edgePair = lath.ledgePair;
    end

    nEdge  = size(edgeAxes,1);
    tiltE  = zeros(nEdge,1);
    coincE = zeros(nEdge,1);
    crestE = nan(nEdge,3);
    azE    = nan(nEdge,1);
    for ii = 1:nEdge
        d  = edgeAxes(ii,:);
        if d(3) < 0, d = -d; end                      % a line: sign irrelevant
        tiltE(ii)  = asind(min(abs(d(3)),1));
        coincE(ii) = cosd(tiltE(ii)) ^ lath.coincExp;
        h = d - d(3)*Z;                               % observable crest = in-surface projection
        if norm(h) > tol
            crestE(ii,:) = h / norm(h);
            azE(ii)      = atan2d(crestE(ii,2), crestE(ii,1));
        end
        edgeAxes(ii,:) = d;
    end

    keepE = find(~isnan(azE));
    if lath.ledgeSelect == "best" && ~isempty(keepE)
        [~,bestIdx] = max(coincE(keepE) .* edgeW(keepE));
        if lath.lathHabit == "rectangular"
            % keep the whole lath, i.e. BOTH rim directions of that rectangle
            keepE = keepE(edgePair(keepE) == edgePair(keepE(bestIdx)));
        else
            keepE = keepE(bestIdx);
        end
    end

    % BEVEL: each rim is a convex edge between the basal broad face (outward
    % normal c_S) and a SIDE face.  The side face contains the rim line and
    % the c-axis, so its normal s is perpendicular to BOTH the rim line and
    % c_S - the lath is a box, so basal and side meet at a right angle and the
    % rounding turns through arcWedge = 90 deg exactly.  The arc therefore
    % runs from c_S round to s and is centred on their BISECTOR, i.e. at
    % psi(c_S) -/+ 45 deg, not on c_S itself.
    %
    % Every rim direction has TWO parallel rims (opposite sides of the
    % rectangle), with side normals +s and -s, so both bisectors occur and the
    % two arcs together cover roughly +/- 90 deg about psi(c_S).  That is why
    % the parallel-to-c stripe comes out with TWO comparable lobes rather than
    % one - as the measured DRP of Fig 2.12 point 4 shows.  rimBias makes the
    % two rims unequal (beta is more etch-resistant on one side than the
    % other); the chapter does not quantify it, so it defaults to 0.
    for ii = keepE(:).'
        ampE = lath.ampBasalEdge * wv * edgeW(ii) * exposeBasal * coincE(ii);
        if ampE <= 0, continue; end
        uE = local_arcBasis(crestE(ii,:));

        % side-face normal: perpendicular to the rim line and to c
        sFace = cross(edgeAxes(ii,:), cS);
        if norm(sFace) < tol, continue; end
        sFace = sFace / norm(sFace);

        if edgeKind(ii) == 2
            kindStr = "lil";     % side face along the lattice invariant line
        else
            kindStr = "prism";   % prismatic side face
        end

        for rim = [1 -1]
            wRim = 1 + rim*lath.rimBias;
            if wRim <= 0, continue; end
            sr  = rim * sFace;
            bis = cS + sr;
            if norm(bis) < tol, continue; end
            bis = bis / norm(bis);

            f = local_blank();
            f.name     = "basalEdge";
            f.kind     = "ridge";
            f.variant  = v;
            f.axis     = edgeAxes(ii,:);
            f.crest    = crestE(ii,:);
            f.crestAz  = azE(ii);
            f.stripeAz = mod(azE(ii) + 90, 180);
            f.tilt     = tiltE(ii);
            f.coinc    = coincE(ii);
            f.expose   = exposeBasal;
            f.amp      = ampE * wRim;
            f.arcU     = uE;
            f.arcCentre   = local_psi(bis, uE);
            f.arcCentreTh = local_arcTheta(f.arcCentre);
            f.arcCentreAz = local_arcAz(f.arcCentre, f.arcU);
            f.arcWedge    = acosd(max(min(dot(cS,sr),1),-1));   % 90 deg
            f.arcSpan     = f.arcWedge + lath.bevelSpan;
            f.edgeKind    = kindStr;
            f.rim         = rim;
            f.normal      = sr;                                 % the side face
            f.note = "basal/" + kindStr + " rim, axis PERP to c -> stripe || c";
            k = k + 1;  R(k,1) = f;
        end
    end

    % =====================================================================
    % FEATURE 3  "basalPlateau":  the flat basal terrace as a plane mirror
    % ---------------------------------------------------------------------
    % A DRP cell (theta,phi) is bright when the microfacet normal that mirrors
    % the illumination into the fixed camera equals the facet normal.  That
    % probe normal is thph2vec(45 + theta/2, phi), i.e. it sits at elevation
    % 45 + theta/2.  The basal normal is at elevation 90 - Phi, so the basal
    % specular point appears at
    %       theta = 2*((90 - Phi) - 45) = 90 - 2*Phi
    % which is exactly the "|90 - 2*Phi|, or 14.4 deg for point 4" quoted in
    % section 2.3.2.  Negative theta = the mirror can never face the camera.
    % Shallow theta is damped because neighbouring peaks shadow the terrace.
    % =====================================================================
    nB     = cS;
    thB    = 2*(asind(min(max(nB(3),-1),1)) - 45);
    shadB  = local_shadow(thB, lath.shadowTh);
    ampB   = lath.ampBasalPlateau * wv * exposeBasal * shadB;
    if ampB > 0
        f = local_blank();
        f.name    = "basalPlateau";
        f.kind    = "peak";
        f.variant = v;
        f.normal  = nB;
        f.peakTh  = thB;
        f.peakAz  = cAz;
        f.expose  = exposeBasal;
        f.shadow  = shadB;
        f.amp     = ampB;
        f.note    = "basal terrace plane mirror -> point maximum at theta = 90 - 2*Phi";
        k = k + 1;  R(k,1) = f;
    end

    % =====================================================================
    % FEATURE 4  "prismPlateau":  the prism ledge terrace as a plane mirror
    % ---------------------------------------------------------------------
    % The R2 intensity marked for point 1 in Fig 2.15.  Same specular rule as
    % feature 3, with the prism normal in place of the basal normal:
    %       theta = 2*(elevation of the prism normal - 45)
    % Amplitude is the projected terrace area (n.Z) plus the pyramidal term:
    % as the terrace RIDGE loses coincidence with the surface (coincT -> 0)
    % the semi-cylinder becomes a semi-spherical peak, so stripe intensity is
    % traded for a point maximum.  pyramidFrac sets how much is traded.
    %
    % prismSelect = "all" (default) offers every upward-facing prism plane as
    % a candidate terrace and lets the geometry select: a peak only appears if
    % its theta = 2*(elevation - 45) falls inside the measured window.  Worked
    % example, point 1 of Fig 2.12 (c almost in the surface, so the prism
    % normals lie on a near-vertical great circle): the near-vertical normal
    % gives theta = 90 deg (unreachable, it IS the macroscopic surface), while
    % the normals 30 deg from Z give theta = 30 deg - reachable, and their
    % azimuth is cAz +/- 90, i.e. the SAME diameter as the terrace stripe.
    % That is precisely how R1.2 and R2 are drawn on top of each other in the
    % Fig 2.15 DRP for point 1.
    % =====================================================================
    keepP = 1:size(prismS,1);
    if lath.prismSelect == "best"
        [~,bestIdx] = max(prismS(:,3) .* lath.prismW);
        keepP = bestIdx;
    end

    for ii = keepP
        nP    = prismS(ii,:);
        thP   = 2*(asind(min(max(nP(3),-1),1)) - 45);
        shadP = local_shadow(thP, lath.shadowTh);
        ampP  = (lath.ampPrismPlateau + lath.pyramidFrac*(1 - coincT)) * ...
                wv * lath.prismW(ii) * (max(nP(3),0)^lath.exposeExp) * shadP;
        if ampP <= 0, continue; end
        f = local_blank();
        f.name    = "prismPlateau";
        f.kind    = "peak";
        f.variant = v;
        f.normal  = nP;
        f.peakTh  = thP;
        f.peakAz  = atan2d(nP(2), nP(1));
        f.expose  = max(nP(3),0)^lath.exposeExp;
        f.shadow  = shadP;
        f.coinc   = coincT;
        f.amp     = ampP;
        f.note    = "prism ledge terrace mirror / pyramidal apex -> point maximum";
        k = k + 1;  R(k,1) = f;
    end
end
end

% =========================================================================
function u = local_arcBasis(t)
% In-surface basis vector of a stripe's great circle.  The circle whose pole
% is the horizontal crest t is spanned by {Z, u} with u = Z x t, so u is
% horizontal and its azimuth IS the stripe azimuth.  Positions along the
% stripe can then be measured as an angle in the {Z,u} plane.
u = cross([0 0 1], t);
n = norm(u);
if n < 1e-9
    u = [NaN NaN NaN];
else
    u = u / n;
end
end

% =========================================================================
function psi = local_psi(n, u)
% Where a direction n sits along the stripe, as an angle from the DRP centre
% (+Z) measured toward the stripe azimuth.  psi > 0 is the half of the stripe
% at stripeAz, psi < 0 the half at stripeAz + 180.
if any(isnan(u)) || any(isnan(n))
    psi = NaN;
else
    psi = atan2d(n(1)*u(1) + n(2)*u(2) + n(3)*u(3), n(3));
end
end

% =========================================================================
function th = local_arcTheta(psi)
% DRP theta corresponding to a position psi along a stripe.  A cell (theta,phi)
% probes the direction at elevation 45 + theta/2, and psi = 90 - elevation, so
% theta = 2*(45 - |psi|).  Values outside [th_min, th_max] mean the arc centre
% falls outside the measured annulus.
th = 2*(45 - abs(psi));
end

% =========================================================================
function az = local_arcAz(psi, u)
% Which of the two stripe lobes the arc centre lies on, as a FULL 0-360
% azimuth.  It has to be derived from the arc basis u, not from the .stripeAz
% field: a stripe is a line, so stripeAz is deliberately folded onto [0,180),
% which throws away exactly the half-line information needed here.  psi > 0 is
% the lobe at azimuth(u), psi < 0 the lobe 180 deg from it.
if isnan(psi) || any(isnan(u))
    az = NaN;
    return
end
az = atan2d(u(2), u(1));
if psi < 0
    az = az + 180;
end
az = mod(az, 360);
end

% =========================================================================
function s = local_shadow(th, shadowTh)
% Specular points are only reachable at positive DRP theta, and shallow
% elevations are progressively lost to shadowing by neighbouring peaks
% (section 2.3.2, point 4).  shadowTh = 0 disables the damping.
if th <= 0
    s = 0;
elseif shadowTh > 0
    s = min(1, th / shadowTh);
else
    s = 1;
end
end

% =========================================================================
function f = local_blank()
% One prototype struct so every element of R has the same fields in the same
% order (required for struct-array assignment).
f = struct( ...
    'name',     "", ...
    'kind',     "", ...
    'variant',  1, ...
    'axis',     [NaN NaN NaN], ...
    'crest',    [NaN NaN NaN], ...
    'crestAz',  NaN, ...
    'stripeAz', NaN, ...
    'tilt',     NaN, ...
    'coinc',    NaN, ...
    'arcU',       [NaN NaN NaN], ...
    'arcCentre',  NaN, ...
    'arcCentreTh',NaN, ...
    'arcCentreAz',NaN, ...
    'arcWedge',   NaN, ...
    'arcSpan',    NaN, ...
    'edgeKind',   "", ...
    'rim',        NaN, ...
    'normal',   [NaN NaN NaN], ...
    'peakTh',   NaN, ...
    'peakAz',   NaN, ...
    'expose',   NaN, ...
    'shadow',   NaN, ...
    'amp',      0, ...
    'note',     "");
end
