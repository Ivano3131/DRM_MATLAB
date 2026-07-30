function exp_para = lathGeometryTi64(exp_para, options)
% lathGeometryTi64  Set up the wrought-Ti64 alpha-lath etching model.
%
%   exp_para = lathGeometryTi64(exp_para)
%   exp_para = lathGeometryTi64(exp_para, name=value, ...)
%
% PURPOSE ------------------------------------------------------------------
% Builds exp_para.lath: the numeric, MTEX-free description of the four
% reflecting features that Kroll's-etched wrought Ti64 presents to the DRM,
% following section 2.3.2 of the thesis chapter ("Directional reflectance
% patterns", Figs 2.12 - 2.15).  lathReflectors / DRPsim_Ti64 consume it.
%
% THE PHYSICAL MECHANISM BEING ENCODED ------------------------------------
% Wrought Ti64 is a stack of ~2 um alpha laths separated by a ~150 nm fine
% beta matrix (Fig 2.14a).  Kroll's reagent attacks alpha much faster than
% beta, and within alpha the BASAL plane etches slowest of all.  The etched
% surface is therefore a lamellar corrugation, and the crests of that
% corrugation are ROUNDED (semi-cylindrical), because the edges between
% differently-etching planes curve (Fig 2.15).
%
% A semi-cylindrical crest along a horizontal axis t has surface normals
% filling the whole great circle perpendicular to t.  It therefore mirrors
% the fixed overhead camera along that entire great circle -> a STRIPE in
% the DRP.  Because t lies in the sample surface, the great circle contains
% +Z and the stripe is a DIAMETER through the DRP centre (observed: "the
% stripes consistently pass through the DRP centroid").
%
%   STRIPE AZIMUTH = (azimuth of the crest t) + 90 deg.
%
% Two crest families exist, at right angles to each other:
%
%   (1) "terraceRidge"  - the ridge that connects the two exposed PRISM
%       facets of a lath, e.g. (10-10) and (11-20) at point 1 of Fig 2.12.
%       Every prism plane contains c, so any prism-prism intersection line
%       is PARALLEL TO THE c-AXIS.  Crest || c  ->  DRP stripe PERPENDICULAR
%       to the c-axis trace.  Needs prism planes to be exposed, i.e. large
%       Phi (c near the surface).  This is the dominant feature at point 1
%       (Phi = 87.9 deg).
%
%   (2) "basalEdge"      - the rounded edge between the chemically resistant
%       beta matrix and the slow-etching BASAL terrace, e.g. point 4 of
%       Fig 2.12.  That edge lies in the basal plane, so it is PERPENDICULAR
%       TO THE c-AXIS  ->  DRP stripe PARALLEL to the c-axis trace.  Needs
%       basal to be exposed, i.e. small Phi.  Dominant at point 4
%       (Phi = 37.8 deg).
%
% Both crests are always present; only their relative amplitude changes with
% Phi, which is exactly the observation that a DRP can show two stripes at
% different intensity magnitudes (section 2.3.2, "lower intensity stripe
% features", 57.6% multi-stripe pixels).
%
% Two point-like (plane-mirror) features complete the model:
%
%   (3) "basalPlateau"  - the flat basal terrace behaves as a plane mirror
%       and gives a POINT maximum, reachable only at illumination elevation
%       theta = 90 - 2*Phi (14.4 deg for point 4, as quoted in the text).
%   (4) "prismPlateau"  - same for the exposed prism ledge terrace (the R2
%       intensity marked for point 1 in Fig 2.15).  Its weight also absorbs
%       the "pyramidal" case: when the terrace ridge line does not coincide
%       with the sample surface it degenerates from a semi-cylinder into a
%       semi-spherical peak, trading stripe intensity for a point maximum
%       (the phi2 effect described at the end of section 2.3.2).
%
% CRYSTALLOGRAPHY ---------------------------------------------------------
% All crystal-frame vectors are produced ONCE here with MTEX (correct
% hexagonal metric and 6/mmm symmetry, same crystalSymmetry object as the
% orientation grid and the IPF colouring), and stored as plain numeric
% arrays so that the per-orientation code (lathReflectors) and the parfor
% dictionary loop need no MTEX object.
%
% OPTIONS -----------------------------------------------------------------
% Crystallography
%   lathHabit      shape of the alpha lath, which fixes how many distinct
%                  edges the basal broad face has:
%                  "rectangular" (DEFAULT, and what wrought Ti64 actually
%                     forms) - the lath is a 6-faced platelet: 2 basal broad
%                     faces, 2 PRISMATIC side faces, and 2 side faces
%                     containing the LATTICE INVARIANT LINE (the lath growth
%                     direction, Fig 2.14).  The basal face is therefore an
%                     approximate rectangle with just TWO edge directions -
%                     one prismatic, one along the lattice invariant line.
%                  "hexagonal" - the old behaviour: every symmetry equivalent
%                     of ledgeUVTW is offered, which is a hexagonal platelet
%                     (6 prismatic side faces + 2 basal = 8 faces).  Kept for
%                     comparison; it is NOT the wrought Ti64 habit.
%   lilUVW         1 x 3 <uvw> of the lattice invariant line in the PARENT
%                  BETA frame (default [-3 3 5], from Fig 2.14b).
%                  IMPORTANT: it is interpreted in the beta variant of that
%                  same figure, [0001]alpha // (110)beta and
%                  <11-20>alpha // [-111]beta, which is NOT variant 1 of
%                  Table 2.1.  The indices must lie in that (110)beta plane
%                  or the line will not come out inside the basal plane; the
%                  code checks this and errors if it does not.
%   ledgeUVTW      N x 4 <uvtw> candidate BASAL-PLANE directions produced by
%                  the PRISMATIC side faces.  Default <10-10>, the ledge
%                  terrace direction named for point 4.  Add [1 1 -2 0] to
%                  also allow <11-20> edges.
%   ledgeWeights   N x 1 relative weight per ledge family (default 1).
%   edgeWeights    1 x 2 [prismatic rim, lattice-invariant-line rim] relative
%                  weight (default [0.3 1]).  This is the lath ASPECT RATIO:
%                  a lath grows long along the lattice invariant line and is
%                  narrow across it (Fig 2.14a), so the two rims running along
%                  the LIL are much longer than the two prismatic end rims and
%                  contribute proportionally more reflecting length.  Use
%                  [1 1] for an equiaxed platelet.  Note the two rim
%                  directions are 75.6 deg apart, so a single lath can produce
%                  a NON-ORTHOGONAL pair of stripes on its own - an
%                  alternative to the variant-superposition explanation.
%   rimBias        -1..1, relative weight of the two PARALLEL rims of a basal
%                  face (default 0 = equal).  The two rims are not physically
%                  equivalent - one meets the beta lamella above, the other
%                  below, and beta is the more etch-resistant phase - but the
%                  chapter does not quantify it, so it is left as a fitted
%                  number.  Non-zero values make the parallel-to-c stripe
%                  one-sided.
%   ledgeSelect    "best"  (default) only the single ledge direction that
%                          best coincides with the sample surface is
%                          realised -> ONE stripe, phi2-sensitive.  This
%                          models lath variant selection.
%                  "trace" use the exact basal surface trace c x Z -> ONE
%                          stripe exactly parallel to the c-axis trace,
%                          phi2-independent.
%                  "all"   every symmetry equivalent -> an asterisk of
%                          stripes (diagnostic only).
%   prismHKIL      P x 4 {hkil} prism planes forming the ledge terraces
%                  (feature 4).  Default {10-10} and {11-20}, the pair
%                  named in Fig 2.15.
%   prismWeights   P x 1 relative weight per prism family (default 1).
%   prismSelect    "all" (default) every upward-facing prism plane is a
%                  candidate terrace; only those whose specular direction
%                  falls inside the measured theta window actually show up,
%                  so the geometry does the selecting.  "best" keeps only the
%                  most upward-facing one (usually the macroscopic surface
%                  itself, whose specular point is unreachable - use with
%                  care).
%
% Amplitudes (relative; the DRP is max-normalised at the end)
%   ampTerraceRidge  stripe PERPENDICULAR to c   (default 1.00)
%   ampBasalEdge     stripe PARALLEL to c        (default 1.00)
%   ampBasalPlateau  basal specular point        (default 0.25)
%   ampPrismPlateau  prism specular point        (default 0.25)
%   pyramidFrac      extra point-peak amplitude gained as the terrace ridge
%                    loses coincidence with the surface (default 0.50)
%
% Shape / weighting exponents
%   exposeExp    exponent on the projected-area factors cos(Phi) / sin(Phi)
%                that say how much basal / prism area the surface exposes
%                (default 1).
%   coincExp     exponent on cos(tilt), the factor that says how well a
%                crest LINE lies in the sample surface: tilt = 0 -> a full
%                semi-cylinder (stripe), tilt -> 90 deg -> a pyramidal peak
%                (default 2).  Raise it to sharpen the phi2 sensitivity.
%   stripeWidth  cross-stripe half-width, deg (default 6).
%   peakWidth    specular peak half-width, deg (default 7).
%   bevelSpan    EXTRA arc beyond the facet wedge, deg (default 0).
%                A real etched crest is a BEVEL, not a full semi-cylinder.
%                Rounding a corner sweeps the surface normal from one face
%                round to the other and NO FURTHER - the flats bound it - so
%                the lit arc is the wedge angle itself, and the geometry fixes
%                it with no free parameter.  Per feature lathReflectors sets
%                    arcSpan   = arcWedge + bevelSpan
%                    arcCentre = bisector of the two facets it joins
%                with arcWedge = 30 deg for the prism-prism terrace ridge and
%                exactly 90 deg for a basal rim (the lath is a box, so basal
%                meets a side face at a right angle).  Leave bevelSpan at 0
%                unless a measured stripe is clearly LONGER than the facet
%                wedge allows.  Use bevelSpan >= 360 to switch the gate off
%                entirely and recover the idealised symmetric stripe.
%   bevelBias    signed tilt added to the TERRACE RIDGE arc centre, deg
%                (default 0).  The two prism facets that crest joins normally
%                straddle the surface normal, so their bisector sits close to
%                +Z and the model predicts a nearly SYMMETRIC perpendicular-
%                to-c stripe.  The real crest is not symmetric - one flank is
%                the etch-resistant beta lamella (Fig 2.15 shows an asymmetric
%                sawtooth; Fig 2.14b puts the lath broad face ~14 deg off the
%                ledge terrace normal) - but that flank cannot be derived from
%                the alpha orientation alone without committing to a parent
%                beta variant.  So it is left here as a fitted number: expect
%                order 10-20 deg, and the SIGN chooses which lobe brightens.
%                The basal edge does not need this; its arc centre is pinned
%                by the basal terrace normal at psi = +/- Phi.
%   bevelSoft    soft shoulder of that arc, deg (default 30).  The rounding is
%                not sharply bounded, so intensity fades over this angle
%                rather than being cut off; 0 gives a hard cut.
%                NOTE the DRP only samples half-angles 45 + theta/2, i.e.
%                |psi| in [12.5, 40] deg for a 10-65 deg theta range, so the
%                two lobes of a stripe are separated by 25-80 deg of arc.
%                Spans much wider than that cannot discriminate them and the
%                stripe stays symmetric.
%   stripeTaper  >=0 exponent on (n.Z) along the stripe (default 2).  This is
%                the microfacet foreshortening pair: the crest element with
%                normal n receives illumination in proportion to n.L and shows
%                the camera an area in proportion to n.Z, and for a specular
%                half-vector n.L = n.Z - hence exponent 2.  It brightens the
%                stripe toward high theta, as observed.  0 = off.
%   shadowTh     specular points below this DRP theta are damped by
%                shadowing from neighbouring peaks, deg (default 20;
%                0 disables).  Section 2.3.2 invokes this for point 4.
%   combine      "max" (default, brightest feature wins - same rule as
%                DRPsim/DRPsim_hcp) | "sum" (incoherent addition).
%
% Burgers-OR variant superposition (EXPERIMENTAL, default off)
%   variantSuperpose  true adds the DRP of the sibling alpha variants that
%                     grow from the same parent beta grain, which section
%                     2.3.2 blames for NON-ORTHOGONAL stripes inside one
%                     DRM pixel.  See borSiblingRotations for the caveats.
%   variantWeight     relative weight of the siblings (default 0.35).
%   variantCount      how many of the 12 BOR variants to include (default 12).
%
%   selfCheck    true prints a numerical check that the hand-rolled
%                EulerRotate agrees with MTEX orientation.byEuler for this
%                crystalSymmetry (guards the convention risk).
%
% OUTPUT ------------------------------------------------------------------
%   exp_para.lath  struct with the numeric crystal-frame vectors
%                    .cAxis        1 x 3  c = [0001] direction
%                    .basalNormal  1 x 3  (0001) plane normal (== cAxis)
%                    .ledgeDirs    N x 3  candidate basal-plane edge lines
%                    .ledgeW       N x 1
%                    .prismNormals P x 3  prism plane normals
%                    .prismW       P x 1
%                    .variantRot   cell of 3 x 3 crystal-frame rotations
%                    .variantW     matching weights
%                  plus every model parameter listed above.
%
% See also lathReflectors, DRPsim_Ti64, checkTi64Mechanism, makeDRPdic_Ti64.
% -------------------------------------------------------------------------
arguments
    exp_para struct
    % --- crystallography -------------------------------------------------
    options.lathHabit    (1,1) string = "rectangular"
    options.lilUVW       (1,3) double = [-3 3 5]
    options.edgeWeights  (1,2) double = [0.3 1]
    options.rimBias      (1,1) double = 0
    options.ledgeUVTW    double = [1 0 -1 0]
    options.ledgeWeights double = []
    options.ledgeSelect  (1,1) string = "best"
    options.prismHKIL    double = [1 0 -1 0; 1 1 -2 0]
    options.prismWeights double = []
    options.prismSelect  (1,1) string = "all"
    % --- amplitudes ------------------------------------------------------
    options.ampTerraceRidge (1,1) double = 1.00
    options.ampBasalEdge    (1,1) double = 1.00
    options.ampBasalPlateau (1,1) double = 0.25
    options.ampPrismPlateau (1,1) double = 0.25
    options.pyramidFrac     (1,1) double = 0.50
    % --- shape / weighting ----------------------------------------------
    options.exposeExp   (1,1) double = 1
    options.coincExp    (1,1) double = 2
    options.stripeWidth (1,1) double = 6
    options.peakWidth   (1,1) double = 7
    options.bevelSpan   (1,1) double = 0
    options.bevelBias   (1,1) double = 0
    options.bevelSoft   (1,1) double = 30
    options.stripeTaper (1,1) double = 2
    options.shadowTh    (1,1) double = 20
    options.combine     (1,1) string = "max"
    % --- Burgers-OR sibling variants (experimental) ----------------------
    options.variantSuperpose (1,1) logical = false
    options.variantWeight    (1,1) double  = 0.35
    options.variantCount     (1,1) double  = 12
    options.csBeta = []
    options.selfCheck        (1,1) logical = false
end

% ---- crystal symmetry ---------------------------------------------------
% Canonical alpha-Ti: a = 2.95 A, c = 4.68 A -> c/a = 1.586.  The SAME object
% must drive the orientation grid and the IPF colouring, so it is taken from
% exp_para when present.
if ~isfield(exp_para,'crystal') || isempty(exp_para.crystal)
    exp_para.crystal = crystalSymmetry('6/mmm',[2.95 2.95 4.68], ...
        'X||a*','Y||b','mineral','Ti-alpha');
end
cs = exp_para.crystal;

% parent beta, needed to carry the lattice invariant line across the Burgers OR
if isempty(options.csBeta)
    csBeta = crystalSymmetry('m-3m',[3.31 3.31 3.31],'mineral','Ti-beta (bcc)');
else
    csBeta = options.csBeta;
end

if ~ismember(options.ledgeSelect,["best","trace","all"])
    error('lathGeometryTi64:ledgeSelect','ledgeSelect must be "best", "trace" or "all".');
end
if ~ismember(options.prismSelect,["best","all"])
    error('lathGeometryTi64:prismSelect','prismSelect must be "best" or "all".');
end
if ~ismember(options.combine,["max","sum"])
    error('lathGeometryTi64:combine','combine must be "max" or "sum".');
end

lath = struct();

% ---- (a) the c-axis / basal normal -------------------------------------
% [0001] as a DIRECTION and (0001) as a PLANE NORMAL must coincide for any
% crystal system; computing both is a free check that the MTEX alignment and
% the hexagonal metric are being applied as expected.
mDir  = Miller(0,0,0,1,cs,'uvtw');   cDir = local_unit(mDir.xyz);
mNrm  = Miller(0,0,0,1,cs,'hkil');   cNrm = local_unit(mNrm.xyz);
if norm(cDir - cNrm) > 1e-6
    warning('lathGeometryTi64:basal', ...
        '[0001] direction and (0001) normal differ by %.2e - check the crystalSymmetry alignment.', ...
        norm(cDir - cNrm));
end
lath.cAxis       = cDir;
lath.basalNormal = cNrm;

% ---- (b) candidate basal/beta EDGE directions (feature 2) --------------
% The rim of the basal terrace lies IN the basal plane, hence perpendicular
% to c.  WHICH directions are available is set by the lath habit:
%
%   an alpha lath in wrought Ti64 is a rectangular PLATELET, not a hexagonal
%   one.  It has 6 faces - 2 basal broad faces, 2 prismatic side faces, and
%   2 side faces containing the lattice invariant line (the growth direction,
%   Fig 2.14) - so its basal face is an approximate rectangle with exactly
%   TWO edge directions.  Offering all six hexagonal <10-10> equivalents (the
%   "hexagonal" habit) would give the lath 8 faces, which is wrong.
[prismEdges, prismEdgeW] = local_symmetriseSet(options.ledgeUVTW, options.ledgeWeights, cs, 'uvtw');
prismEdges = local_keepBasal(prismEdges, lath.cAxis, 'ledgeUVTW');
if isempty(prismEdges)
    error('lathGeometryTi64:noLedge', ...
        'no basal-plane ledge directions left; ledgeUVTW must be <uvtw> with w = 0.');
end
prismEdgeW = prismEdgeW(1:size(prismEdges,1));

switch options.lathHabit
    case "hexagonal"
        % every prismatic equivalent, no lattice invariant line (old model)
        lath.ledgeDirs = prismEdges;
        lath.ledgeW    = options.edgeWeights(1) * prismEdgeW;
        lath.ledgeKind = ones(size(prismEdges,1),1);       % 1 = prismatic
        lath.ledgePair = (1:size(prismEdges,1)).';         % each on its own
        lath.lilAngle       = NaN;
        lath.edgePairAngle  = NaN;

    case "rectangular"
        [lilDirs, lilAngle] = local_lilDirections(cs, csBeta, options.lilUVW, lath.cAxis);
        % Pair each lattice-invariant-line edge with the prismatic edge that
        % is most nearly perpendicular to it: those two directions are the
        % two sides of the (approximate) rectangle, so together they define
        % one lath.  The included angle is reported because it is NOT exactly
        % 90 deg - the LIL is irrational, ~14.4 deg off <11-20>.
        nL = size(lilDirs,1);
        dirs = zeros(0,3); w = zeros(0,1); kind = zeros(0,1); pair = zeros(0,1);
        incl = zeros(nL,1);
        for k = 1:nL
            d = abs(prismEdges * lilDirs(k,:).');
            [dmin,ip] = min(d);
            incl(k) = acosd(min(dmin,1));
            dirs = [dirs; prismEdges(ip,:); lilDirs(k,:)];              %#ok<AGROW>
            w    = [w; options.edgeWeights(1)*prismEdgeW(ip); ...
                       options.edgeWeights(2)];                          %#ok<AGROW>
            kind = [kind; 1; 2];                                         %#ok<AGROW>
            pair = [pair; k; k];                                         %#ok<AGROW>
        end
        lath.ledgeDirs = dirs;
        lath.ledgeW    = w;
        lath.ledgeKind = kind;      % 1 = prismatic side, 2 = LIL side
        lath.ledgePair = pair;      % rows sharing a pair id are one lath
        lath.lilAngle      = lilAngle;
        lath.edgePairAngle = mean(incl);

    otherwise
        error('lathGeometryTi64:lathHabit', ...
            'lathHabit must be "rectangular" or "hexagonal".');
end

% ---- (c) prism planes forming the ledge terraces (feature 4) -----------
[prismNormals, prismW] = local_symmetriseSet(options.prismHKIL, options.prismWeights, cs, 'hkl');
lath.prismNormals = prismNormals;
lath.prismW       = prismW;

% ---- (d) Burgers-OR sibling variants (optional) ------------------------
if options.variantSuperpose
    [Rlist, vlabels, vangles] = borSiblingRotations(cs, count=options.variantCount);
    lath.variantRot   = Rlist;
    lath.variantW     = [1; options.variantWeight*ones(numel(Rlist)-1,1)];
    lath.variantLabel = vlabels;
    lath.variantAngle = vangles;
else
    lath.variantRot   = {eye(3)};
    lath.variantW     = 1;
    lath.variantLabel = "parent alpha variant only";
    lath.variantAngle = 0;
end

% ---- (e) model parameters (copied verbatim so they travel into parfor) --
lath.lathHabit       = options.lathHabit;
lath.rimBias         = options.rimBias;
lath.ledgeSelect     = options.ledgeSelect;
lath.prismSelect     = options.prismSelect;
lath.ampTerraceRidge = options.ampTerraceRidge;
lath.ampBasalEdge    = options.ampBasalEdge;
lath.ampBasalPlateau = options.ampBasalPlateau;
lath.ampPrismPlateau = options.ampPrismPlateau;
lath.pyramidFrac     = options.pyramidFrac;
lath.exposeExp       = options.exposeExp;
lath.coincExp        = options.coincExp;
lath.stripeWidth     = options.stripeWidth;
lath.peakWidth       = options.peakWidth;
lath.bevelSpan       = options.bevelSpan;
lath.bevelBias       = options.bevelBias;
lath.bevelSoft       = options.bevelSoft;
lath.stripeTaper     = options.stripeTaper;
lath.shadowTh        = options.shadowTh;
lath.combine         = options.combine;

exp_para.lath = lath;

% ---- (f) optional convention self-check --------------------------------
% EulerRotate is a hand-rolled active Bunge rotation (specimen =
% Rz(phi1)*Rx(Phi)*Rz(phi2) * crystal).  Everything in this model rests on
% it agreeing with MTEX, which supplies the crystal-frame vectors and the
% IPF colours.  This compares the two on random orientations.
if options.selfCheck
    nTest = 20;
    dev   = zeros(nTest,1);
    for ii = 1:nTest
        eu = [360*rand, 180*rand, 360*rand];
        vHand = local_unit(EulerRotate(lath.cAxis, eu(1), eu(2), eu(3)));
        o     = orientation.byEuler(eu(1)*degree, eu(2)*degree, eu(3)*degree, cs);
        vm    = o * Miller(0,0,0,1,cs,'uvtw');
        vMtex = local_unit([vm.x vm.y vm.z]);
        dev(ii) = min(norm(vHand - vMtex), norm(vHand + vMtex));   % c is an axis
    end
    fprintf('lathGeometryTi64 selfCheck: EulerRotate vs MTEX c-axis, max |dev| = %.2e\n', max(dev));
    if max(dev) > 1e-6
        warning('lathGeometryTi64:convention', ...
            ['EulerRotate and MTEX disagree (%.2e). Every predicted stripe azimuth ', ...
             'will be wrong. Fix EulerRotate or rotate with MTEX instead.'], max(dev));
    end
end

fprintf(['Ti64 lath model ready: %s habit, %d rim direction(s) [%s], ', ...
         '%d prism normal(s) [%s], %d alpha variant(s).\n'], ...
    lath.lathHabit, size(lath.ledgeDirs,1), lath.ledgeSelect, ...
    size(lath.prismNormals,1), lath.prismSelect, numel(lath.variantRot));
if lath.lathHabit == "rectangular"
    fprintf(['  lath = 6 faces (2 basal + 2 prismatic + 2 lattice-invariant-line); ', ...
             'LIL is %.2f deg from <11-20> (Fig 2.14b labels 14.4),\n', ...
             '  and the two basal rim directions meet at %.2f deg ', ...
             '("essentially a rectangle").\n'], lath.lilAngle, lath.edgePairAngle);
end
end

% =========================================================================
function vecs = local_keepBasal(vecs, cAxis, what)
% Keep only directions lying in the basal plane (perpendicular to c), which a
% terrace rim must be.
perp = abs(vecs * cAxis(:));
keep = perp < 1e-6;
if ~all(keep)
    warning('lathGeometryTi64:notBasal', ...
        '%d of %d %s directions are not perpendicular to c and were dropped.', ...
        sum(~keep), numel(keep), what);
end
vecs = vecs(keep,:);
end

% =========================================================================
function [dirs, angToA] = local_lilDirections(csAlpha, csBeta, lilUVW, cAxis)
% The LATTICE INVARIANT LINE, expressed in the alpha crystal frame.
%
% The line is quoted in the parent beta frame (Fig 2.14b gives [-3 3 5]beta),
% so it has to be carried across the Burgers OR.  It must be carried using the
% SAME beta variant the figure uses,
%       (110)beta // (0001)alpha      and     [-111]beta // [11-20]alpha
% which is NOT variant 1 of Table 2.1 ((1-10)beta).  Using the wrong variant
% silently produces a line that does not lie in the basal plane at all - it
% comes out ~50 deg off - so the result is checked below rather than trusted.
%
% For the default indices this returns 6 distinct axes at 14.42 deg from the
% <11-20> directions, reproducing the 14.4 deg angle labelled in Fig 2.14b.
p = Miller(1,1,0,csBeta,'hkl');       % (110)beta  // (0001)alpha
d = Miller(-1,1,1,csBeta,'uvw');      % [-111]beta // [11-20]alpha
o = orientation.map(p, Miller(0,0,0,1,csAlpha,'hkil'), ...
                    d, Miller(1,1,-2,0,csAlpha,'UVTW'));

if abs(dot(lilUVW(:).', [1 1 0])) > 1e-9
    % NB error() requires SCALAR formatted arguments, unlike fprintf
    error('lathGeometryTi64:lilVariant', ...
        ['lilUVW = [%g %g %g] does not lie in the (110)beta plane, so it is ', ...
         'indexed in a different beta variant than the one used to interpret ', ...
         'it. Re-index it into the (110)beta // (0001)alpha variant.'], ...
        lilUVW(1), lilUVW(2), lilUVW(3));
end

L  = o * Miller(lilUVW(1), lilUVW(2), lilUVW(3), csBeta, 'uvw');
v  = [L.x L.y L.z];
v  = v / norm(v);

% must land in the basal plane; anything else means the OR/variant is wrong
offBasal = 90 - acosd(min(abs(dot(v, cAxis(:).')),1));
if abs(offBasal) > 1e-3
    error('lathGeometryTi64:lilNotBasal', ...
        ['the lattice invariant line came out %.2f deg away from the basal ', ...
         'plane; it must lie in it. Check lilUVW and the beta variant.'], offBasal);
end

% all symmetry equivalents, +/- collapsed onto one axis each
mL  = Miller(vector3d(v(1),v(2),v(3)), csAlpha);
ms  = symmetrise(mL);
xyz = local_canonicalSign(normr(ms.xyz));
dirs = uniquetol(xyz,1e-4,'ByRows',true,'DataScale',1);

% report the angle to the nearest a-axis (the Fig 2.14b label, 14.4 deg)
aM  = Miller(1,1,-2,0,csAlpha,'UVTW');
aS  = symmetrise(aM);
aXY = normr(aS.xyz);
angToA = min(acosd(min(abs(aXY * dirs(1,:).'),1)));
end

% =========================================================================
function [vecs, w] = local_symmetriseSet(idx4, wIn, cs, convention)
% Expand N x 4 Miller-Bravais rows into the full symmetry-equivalent set of
% unit vectors in the crystal Cartesian frame.  +v and -v describe the same
% line / the same plane, so they are collapsed onto one representative.
if size(idx4,2) ~= 4
    error('lathGeometryTi64:badIndices','Miller-Bravais input must have 4 columns.');
end
N = size(idx4,1);
if isempty(wIn)
    wIn = ones(N,1);
else
    wIn = wIn(:);
    if numel(wIn) ~= N
        error('lathGeometryTi64:badWeights','need one weight per Miller-Bravais row.');
    end
end

vecs = zeros(0,3);
w    = zeros(0,1);
for k = 1:N
    m   = Miller(idx4(k,1),idx4(k,2),idx4(k,3),idx4(k,4),cs,convention);
    ms  = symmetrise(m);
    xyz = normr(ms.xyz);
    xyz = local_canonicalSign(xyz);
    % 'DataScale',1 is essential: uniquetol otherwise scales the tolerance by
    % the data magnitude COLUMN BY COLUMN, and a column that is numerically
    % zero (e.g. the z component of basal-plane directions, ~1e-16) then gets
    % an effective tolerance of ~1e-22, so round-off copies of the same vector
    % survive as separate "symmetry equivalents".  With the default scaling
    % <10-10> comes out as 12 directions instead of 3.
    xyz = uniquetol(xyz,1e-4,'ByRows',true,'DataScale',1);
    vecs = [vecs; xyz];                        %#ok<AGROW>
    w    = [w;    wIn(k)*ones(size(xyz,1),1)]; %#ok<AGROW>
end
% NOTE: deliberately NOT sum-normalised.  Features are combined with max(),
% so each candidate must carry its own absolute weight; dividing by the
% family size would make a 6-member family six times dimmer than a 1-member
% one for no physical reason.
end

% =========================================================================
function xyz = local_canonicalSign(xyz)
tol = 1e-8;
for r = 1:size(xyz,1)
    v   = xyz(r,:);
    idx = find(abs(v) > tol,1,'first');
    if ~isempty(idx) && v(idx) < 0
        xyz(r,:) = -v;
    end
end
end

% =========================================================================
function v = local_unit(v)
v = double(v(:)).';
n = norm(v);
if n < eps
    error('lathGeometryTi64:zeroVector','zero-length crystal vector.');
end
v = v / n;
end
