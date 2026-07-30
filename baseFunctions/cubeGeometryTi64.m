function exp_para = cubeGeometryTi64(exp_para, options)
% cubeGeometryTi64  The alpha lath as a box with reflecting beveled edges.
%
%   exp_para = cubeGeometryTi64(exp_para)
%   exp_para = cubeGeometryTi64(exp_para, dims=[1 1 1], etchFrac=0.3, ...)
%
% THE MODEL IN ONE PARAGRAPH ----------------------------------------------
% An alpha lath in wrought Ti64 is treated as a BOX with three orthogonal
% face-pairs:
%
%     (0001)    basal, normal e_c   <- etches slowest
%     {10-10}   the lath BROAD face, normal e_b
%     {11-20}   the lath END face,   normal e_p
%
% Those three normals are mutually perpendicular in an HCP lattice (the
% {10-10} and {11-20} normals lie in the basal plane 90 deg apart), so the
% object really is a rectangular box, and with dims = [1 1 1] a cube.
%
% Only the 12 EDGES reflect.  Etching rounds each edge into a bevel, and a
% rounded edge with axis t acts as a convex cylindrical mirror: it fans the
% light out along the great circle perpendicular to t, which is a STRIPE in
% the DRP.  Flat faces are ignored entirely - a polished plane mirror only
% lights one DRP cell, which is not what a DRP of etched Ti64 looks like.
%
% PARTIAL ETCHING is a single cut plane.  After the box is rotated into the
% sample frame, only the top `etchFrac` of its vertical extent is taken to
% stand proud of the etched surface; an edge contributes in proportion to the
% fraction of its length that survives above the cut.  Nothing else selects
% which edges are visible.
%
% That one rule reproduces the section-2.3.2 mechanism without asserting it:
%
%   Phi ~ 88 deg (Fig 2.12 point 1) - c lies in the surface, the box stands on
%     edge, and the exposed top edges are the four PARALLEL to c (the
%     broad-face / end-face ridges).  Their stripes are PERPENDICULAR to c.
%
%   Phi ~ 38 deg (Fig 2.12 point 4) - the basal face is near-horizontal and
%     forms the top, so the exposed edges are its rim, which is PERPENDICULAR
%     to c.  Their stripes are PARALLEL to c.
%
% EDGE FAMILIES -----------------------------------------------------------
% Edges are named by the two faces they join, and grouped by their axis:
%
%   axis 1  || e_c   broadEnd     broad face  ^ end face    -> stripe PERP to c
%   axis 2  || e_b   basalEnd     basal face  ^ end face    -> stripe || to c
%   axis 3  || e_p   basalBroad   basal face  ^ broad face  -> stripe || to c
%
% INPUTS ------------------------------------------------------------------
%   exp_para  struct; needs .crystal (an MTEX crystalSymmetry).  If absent,
%             the alpha-Ti default a = 2.95 A, c = 4.68 A is used.
%
% OPTIONS -----------------------------------------------------------------
%   dims        [Lc Lb Lp] box side lengths along [e_c e_b e_p].  [1 1 1] is
%               a cube; a plate-like lath is e.g. [0.4 1 1].  Only the RATIOS
%               matter.  Longer edges reflect more light, so dims also sets
%               the relative edge weights.
%   etchFrac    0..1, fraction of the rotated box height left standing proud
%               of the etched surface (default 0.30).  1 exposes every edge
%               (symmetric asterisk); small values expose only the topmost
%               family.  THIS is the knob that controls which stripe wins.
%   exposeExp   exponent on the exposed length fraction (default 1).  Raise
%               it to make the etch selection sharper.
%   edgeAmp     intrinsic reflectivity of the three edge families, indexed by
%               the axis the edge runs along (default [1 1 1]):
%                 (1) edges || e_c  = "broadEnd",   broad ^ end   (prism^prism)
%                 (2) edges || e_b  = "basalEnd",   basal ^ end
%                 (3) edges || e_p  = "basalBroad", basal ^ broad
%               so (2) and (3) together are the BASAL RIM.  Raising them above
%               (1) says the slow-etching basal plane leaves a crisper bevel,
%               which is the one chemistry argument the model needs; the rest
%               of the work should be done by the geometry.
%   stripeWidth cross-stripe Cauchy half-width, deg (default 6).
%   bevelSpan   EXTRA arc beyond the 90 deg facet wedge, deg (default 0).
%               Rounding an edge sweeps its normal from one face round to the
%               other and no further, so the lit arc IS the wedge - there is
%               no free parameter unless a measured stripe is longer than the
%               wedge allows.  >= 270 removes the gate (symmetric stripe).
%   bevelSoft   raised-cosine shoulder on that arc, deg (default 30).
%   stripeTaper exponent on (v.Z) along the stripe (default 2 = the
%               microfacet foreshortening pair, illumination x visible area).
%   combine     "max" brightest edge wins (default) | "sum" incoherent
%               addition.
%   selfCheck   true (default) verifies the triad is orthonormal, that e_b
%               really is a {10-10} normal, and that EulerRotate agrees with
%               MTEX.
%
% OUTPUT ------------------------------------------------------------------
%   exp_para.cube  a pure-numeric struct (no MTEX objects, so it travels into
%                  parfor):
%     .eC .eB .eP    1x3 crystal-frame unit vectors
%     .dims          1x3
%     .vert          8x3 vertices, crystal frame
%     .edgeV         12x2 vertex index pairs
%     .edgeAxis      12x1 which axis (1,2,3) the edge runs along
%     .edgeN1/.edgeN2 12x3 outward normals of the two faces the edge joins
%     .edgeLenW      12x1 length weight, normalised to max 1
%     .edgeAmp       12x1 intrinsic reflectivity
%     .family        12x1 string
%     plus every scalar option above.
%
% See also cubeReflectors, DRPsim_cube, checkCubeMechanism, makeDRPdic_cube.
% -------------------------------------------------------------------------
arguments
    exp_para struct
    options.dims (1,3) double = [1 1 1]
    options.etchFrac (1,1) double {mustBeGreaterThan(options.etchFrac,0), ...
                                   mustBeLessThanOrEqual(options.etchFrac,1)} = 0.30
    options.exposeExp (1,1) double = 1
    options.edgeAmp (1,3) double = [1 1 1]
    options.stripeWidth (1,1) double = 6
    options.bevelSpan (1,1) double = 0
    options.bevelSoft (1,1) double = 30
    options.stripeTaper (1,1) double = 2
    options.combine (1,1) string = "max"
    options.selfCheck (1,1) logical = true
end

if ~any(options.combine == ["max" "sum"])
    error('cubeGeometryTi64:combine','combine must be "max" or "sum".');
end
if any(options.dims <= 0)
    error('cubeGeometryTi64:dims','dims must be positive.');
end

if ~isfield(exp_para,'crystal') || isempty(exp_para.crystal)
    exp_para.crystal = crystalSymmetry('6/mmm',[2.95 2.95 4.68], ...
        'X||a*','Y||b','mineral','Ti-alpha');
end
cs = exp_para.crystal;

% ---- (a) the orthonormal crystal-frame triad ----------------------------
% e_c is fixed by the lattice.  e_p is a {11-20} normal, which lies in the
% basal plane.  e_b is then FORCED: it is the only direction perpendicular to
% both, and it comes out on a {10-10} normal - which is why a basal / broad /
% end box can be exactly rectangular in a hexagonal crystal.
eC = local_unit(Miller(0,0,0,1,cs,'hkil').xyz);
eP = local_unit(Miller(1,1,-2,0,cs,'hkl').xyz);
eB = local_unit(cross(eC,eP));
eP = local_unit(cross(eB,eC));      % re-orthogonalise against round-off

cube = struct();
cube.eC = eC;  cube.eB = eB;  cube.eP = eP;
cube.dims = options.dims;

% ---- (b) vertices -------------------------------------------------------
% Basis matrix B has the three axes as ROWS, matching EulerRotate's
% row-vector convention (v_sample = v_crystal * R).
B = [eC; eB; eP];
S = local_signCombos();                       % 8x3 of +-1
cube.vert = 0.5 * (S .* options.dims) * B;    % 8x3, crystal frame

% ---- (c) the 12 edges ---------------------------------------------------
% An edge parallel to axis k joins the two faces whose normals are the OTHER
% two axes; its position is fixed by the signs on those two axes.  Walking
% (k, sm, sn) therefore enumerates all 12 edges exactly once.
famName = ["broadEnd" "basalEnd" "basalBroad"];   % indexed by axis 1,2,3
axes3   = {eC, eB, eP};

edgeV    = zeros(12,2);
edgeAxis = zeros(12,1);
edgeN1   = zeros(12,3);
edgeN2   = zeros(12,3);
edgeLenW = zeros(12,1);
edgeAmp  = zeros(12,1);
family   = strings(12,1);

e = 0;
for k = 1:3
    o = setdiff(1:3,k);          % the two axes that fix this edge's position
    m = o(1); n = o(2);
    for sm = [-1 1]
        for sn = [-1 1]
            e = e + 1;
            s = zeros(1,3);
            s(m) = sm;  s(n) = sn;
            s(k) = -1;  v1 = local_vertexIndex(S,s);
            s(k) = +1;  v2 = local_vertexIndex(S,s);

            edgeV(e,:)    = [v1 v2];
            edgeAxis(e)   = k;
            edgeN1(e,:)   = sm * axes3{m};    % outward normal of face m
            edgeN2(e,:)   = sn * axes3{n};    % outward normal of face n
            edgeLenW(e)   = options.dims(k);
            edgeAmp(e)    = options.edgeAmp(k);
            family(e)     = famName(k);
        end
    end
end

cube.edgeV    = edgeV;
cube.edgeAxis = edgeAxis;
cube.edgeN1   = edgeN1;
cube.edgeN2   = edgeN2;
cube.edgeLenW = edgeLenW / max(edgeLenW);
cube.edgeAmp  = edgeAmp;
cube.family   = family;

% ---- (d) model parameters (plain numerics, so they survive parfor) ------
cube.etchFrac    = options.etchFrac;
cube.exposeExp   = options.exposeExp;
cube.stripeWidth = options.stripeWidth;
cube.bevelSpan   = options.bevelSpan;
cube.bevelSoft   = options.bevelSoft;
cube.stripeTaper = options.stripeTaper;
cube.combine     = options.combine;

exp_para.cube = cube;

% ---- (e) is this box actually indexable? -------------------------------
% If two of the three axes have the SAME side length and the SAME edge
% amplitude, a 90 deg rotation about the third axis maps the box exactly onto
% itself.  That rotation is not in 6/mmm, so it sends a crystal orientation to
% a genuinely DIFFERENT one with an IDENTICAL DRP - the model then cannot tell
% those orientations apart, no matter how good the matcher is.
%
% A literal cube (dims [1 1 1], edgeAmp [1 1 1]) has the full octahedral
% symmetry and is completely unindexable: measured over a 5 deg dictionary,
% every entry has a twin at correlation 1.000 and the median round-trip
% misorientation is 90 deg.  Give the three axes distinct lengths - which a
% real lath has anyway - and the median drops below the grid step.
pairName = ["broad & end (rotation about c)", ...
            "c & end (rotation about the broad-face normal)", ...
            "c & broad (rotation about the end-face normal)"];
pairs = [2 3; 1 3; 1 2];
degen = false(1,3);
for k = 1:3
    p = pairs(k,:);
    degen(k) = abs(diff(options.dims(p))) < 1e-9 && abs(diff(options.edgeAmp(p))) < 1e-9;
end
if any(degen)
    warning('cubeGeometryTi64:degenerate', ...
        ['this box maps onto itself under a 90 deg rotation swapping %s, so ' ...
         'the DRP cannot distinguish those two axes and indexing will be ' ...
         'ambiguous.\nGive the axes distinct dims (a real lath is a plate, ' ...
         'e.g. dims = [1 0.3 1.5]) or distinct edgeAmp.\nRun cubeAmbiguity ' ...
         'to measure how bad it is.'], strjoin(pairName(degen), ' and '));
end

% ---- (f) self-check -----------------------------------------------------
if options.selfCheck
    G   = B * B.';
    dev = max(abs(G - eye(3)),[],'all');
    fprintf('cubeGeometryTi64 selfCheck: triad orthonormality |B*B'' - I| = %.2e\n', dev);
    if dev > 1e-12
        warning('cubeGeometryTi64:triad', ...
            'the basal / broad / end normals are not orthogonal (%.2e).', dev);
    end

    % e_b should sit on a {10-10} normal; if it does not, the box is not the
    % object the model claims it is.
    mB   = symmetrise(Miller(1,0,-1,0,cs,'hkl'));
    xyz  = normr([mB.x(:) mB.y(:) mB.z(:)]);
    angB = min(acosd(min(abs(xyz*eB.'),1)));
    fprintf('cubeGeometryTi64 selfCheck: broad-face normal is %.2e deg from the nearest {10-10}\n', angB);
    if angB > 1e-6
        warning('cubeGeometryTi64:broadFace', ...
            'the broad-face normal is %.3f deg off {10-10}.', angB);
    end

    % EulerRotate is a hand-rolled active Bunge rotation; every predicted
    % stripe azimuth rests on it agreeing with MTEX.
    nTest = 20; devE = zeros(nTest,1);
    for ii = 1:nTest
        eu    = [360*rand, 180*rand, 360*rand];
        vHand = local_unit(EulerRotate(eC, eu(1), eu(2), eu(3)));
        vm    = orientation.byEuler(eu(1)*degree,eu(2)*degree,eu(3)*degree,cs) ...
                * Miller(0,0,0,1,cs,'uvtw');
        vMtex = local_unit([vm.x vm.y vm.z]);
        devE(ii) = min(norm(vHand-vMtex), norm(vHand+vMtex));
    end
    fprintf('cubeGeometryTi64 selfCheck: EulerRotate vs MTEX c-axis, max |dev| = %.2e\n', max(devE));
    if max(devE) > 1e-6
        warning('cubeGeometryTi64:convention', ...
            'EulerRotate and MTEX disagree (%.2e); every stripe azimuth will be wrong.', max(devE));
    end
end

fprintf(['Ti64 cube model ready: dims [%.2f %.2f %.2f] along [c, broad, end], ', ...
         'etchFrac %.2f, 12 beveled edges (4 || c, 8 in the basal plane).\n'], ...
    options.dims, options.etchFrac);
end

% =========================================================================
function u = local_unit(v)
v = double(v(:)).';
u = v / norm(v);
end

% =========================================================================
function S = local_signCombos()
% The 8 vertices, ordered so local_vertexIndex can find one by its signs.
[a,b,c] = ndgrid([-1 1],[-1 1],[-1 1]);
S = [a(:) b(:) c(:)];
end

% =========================================================================
function idx = local_vertexIndex(S, s)
idx = find(all(S == s, 2), 1);
end
