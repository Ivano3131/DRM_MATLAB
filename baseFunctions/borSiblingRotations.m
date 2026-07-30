function [Rlist, labels, angs] = borSiblingRotations(csAlpha, options)
% borSiblingRotations  Crystal-frame rotations to the Burgers-OR sibling
%                      alpha variants of a given alpha grain.
%
%   [Rlist, labels, angs] = borSiblingRotations(csAlpha)
%   [...] = borSiblingRotations(csAlpha, csBeta=..., count=..., verbose=...)
%
% WHY ---------------------------------------------------------------------
% Section 2.3.2 attributes the reported NON-ORTHOGONAL stripes to the Burgers
% Orientation Relationship: up to 12 alpha variants grow from one parent beta
% grain in a closely spaced cluster, so one DRM pixel can superimpose the
% reflectance of several differently-oriented laths.  This function supplies
% the misorientations needed to add those siblings to a simulated DRP.
%
% BOR (equations 2.3 / 2.4 of the chapter)
%       {0001}alpha // {110}beta        <11-20>alpha // <111>beta
% Table 2.1 lists the 12 (plane, direction) pairs; the same list is used here,
% in the same row order, so the printed misorientation angles can be checked
% straight against the table (10.53, 60.00, 60.83, 63.26, 90.00 deg).
%
% CONVENTION --------------------------------------------------------------
% Rlist{i} is the 3 x 3 rotation matrix that takes a vector expressed in
% variant i's crystal frame into variant 1's crystal frame:
%
%       u_variant1 = Rlist{i} * u_varianti          (column vectors)
%
% so a caller that already rotates variant-1 crystal vectors into the
% specimen frame only has to pre-multiply by Rlist{i} to get the sibling's
% vectors.  Rlist{1} is the identity.
%
% LIMITATION (read before trusting the result) ----------------------------
% A single alpha orientation does not determine its parent beta uniquely -
% each alpha variant has SIX possible parents (stated in the chapter).  This
% function returns the sibling set for ONE parent assignment, the one in which
% the indexed grain plays the role of variant 1 of Table 2.1.  A different
% parent assignment would give a different (symmetry-related) sibling set.
% Treat the superposition it enables as exploratory, not as ground truth.
%
% Inputs / options
%   csAlpha   MTEX crystalSymmetry of alpha; its ALIGNMENT matters, because
%             Rlist is expressed in that Cartesian crystal frame.  Pass the
%             same object used everywhere else (exp_para.crystal).
%   csBeta    beta crystalSymmetry (default bcc, a = 3.31 A).
%   count     how many of the 12 variants to return (default 12).
%   verbose   print the misorientation table (default true).
%
% Outputs
%   Rlist   1 x count cell of 3 x 3 rotation matrices (Rlist{1} = eye(3))
%   labels  1 x count string, the Table 2.1 (plane)//(0001) [dir]//[11-20] pair
%   angs    1 x count misorientation angle w.r.t. variant 1, deg
% -------------------------------------------------------------------------
arguments
    csAlpha
    options.csBeta = crystalSymmetry('m-3m',[3.31 3.31 3.31],'mineral','Ti-beta (bcc)')
    options.count (1,1) double = 12
    options.verbose (1,1) logical = true
end

csBeta = options.csBeta;
nWant  = max(1, min(12, round(options.count)));

% ---- Table 2.1, in table order -----------------------------------------
planeBeta = [ 1 -1  0;   1  0 -1;   0  1 -1;      % variants 1-3
              1  1  0;   1  0  1;   0  1 -1;      % variants 4-6
              1  1  0;   1  0 -1;   0  1  1;      % variants 7-9
              1 -1  0;   1  0  1;   0  1  1];     % variants 10-12
dirBeta   = [ 1  1  1;   1  1  1;   1  1  1;
             -1  1  1;  -1  1  1;  -1  1  1;
              1 -1  1;   1 -1  1;   1 -1  1;
              1  1 -1;   1  1 -1;   1  1 -1];

cAlpha = Miller(0,0,0,1,csAlpha,'hkil');     % (0001)
aAlpha = Miller(1,1,-2,0,csAlpha,'UVTW');    % [11-20]

p1 = Miller(planeBeta(1,1),planeBeta(1,2),planeBeta(1,3),csBeta,'hkl');
d1 = Miller(dirBeta(1,1),  dirBeta(1,2),  dirBeta(1,3),  csBeta,'uvw');

% ---- fit the OR once, then generate the variants by beta symmetry -------
% o1 maps beta-frame vectors onto alpha-frame vectors for variant 1.  Variant
% i is variant 1 composed with the cubic operation that carries (p1,d1) onto
% (pI,dI); that operation is a genuine beta symmetry element, so the
% intervariant misorientations come out exact.
o1 = orientation.map(p1, cAlpha, d1, aAlpha);
v  = repmat(o1, 12, 1);
for i = 1:12
    pI   = Miller(planeBeta(i,1),planeBeta(i,2),planeBeta(i,3),csBeta,'hkl');
    dI   = Miller(dirBeta(i,1),  dirBeta(i,2),  dirBeta(i,3),  csBeta,'uvw');
    s_i  = orientation.map(p1, pI, d1, dI);
    v(i) = o1 .* inv(s_i);
end

% ---- crystal-frame rotation from variant i back to variant 1 -----------
% With g_i the specimen orientation of variant i and g_1 that of variant 1,
% g_i = g_1 * Delta_i with Delta_i = v(1) * inv(v(i)); a vector u given in
% variant i's frame has variant-1 coordinates Delta_i * u.
Rlist  = cell(1,nWant);
labels = strings(1,nWant);
angs   = zeros(1,nWant);
% NOTE on the angle: it must be taken from delta itself.  inv(v(1)).*v(i) is
% the other composition order and MTEX types it as a beta-beta misorientation,
% so it gets reduced under CUBIC symmetry and returns 0 or 60 deg for every
% variant - wrong.  angle(delta) reduces under the alpha (hexagonal) symmetry
% and reproduces Table 2.1 exactly, row for row.
for i = 1:nWant
    delta     = v(1) .* inv(v(i));
    Rlist{i}  = matrix(delta);
    angs(i)   = angle(delta) / degree;
    labels(i) = sprintf('v%d: (%d%d%d)//(0001), [%d%d%d]//[11-20]', i, ...
        planeBeta(i,1),planeBeta(i,2),planeBeta(i,3), ...
        dirBeta(i,1),dirBeta(i,2),dirBeta(i,3));
end
Rlist{1} = eye(3);      % exact identity, free of round-off
angs(1)  = 0;

% ---- sanity check against Table 2.1 ------------------------------------
if options.verbose
    fprintf('Burgers-OR sibling variants (misorientation w.r.t. variant 1):\n');
    for i = 1:nWant
        fprintf('  %-46s  %6.2f deg\n', labels(i), angs(i));
    end
    fprintf('  distinct angles: ');
    fprintf('%.2f  ', unique(round(angs(2:end),2)));
    fprintf('\n  Table 2.1 expects: 10.53  60.00  60.83  63.26  90.00\n');
end
end
