%% burgers_variants_fixed.m
%  Misorientation of the 12 Burgers-OR alpha variants w.r.t. variant 1,
%  in the SAME row order as the table.
%
%  Key idea: fit the OR ONCE (variant 1), then build every other variant by
%  composing with the beta symmetry operation s_i that maps variant 1's
%  (plane,direction) onto row i's. Intervariant misorientations are then exact.
% -----------------------------------------------------------------------------
 
%% Crystal symmetries
csBeta  = crystalSymmetry('m-3m', [3.31 3.31 3.31], 'mineral','Ti-beta (bcc)');
csAlpha = crystalSymmetry('6/mmm',[2.95 2.95 4.68], ...
            'X||a','Y||b*','Z||c*', 'mineral','Ti-alpha (hcp)');   % c/a = 1.587
 
%% Table definitions
planeBeta = [ 1 -1  0;   1  0 -1;   0  1 -1;    % 1-3
              1  1  0;   1  0  1;   0  1 -1;    % 4-6
              1  1  0;   1  0 -1;   0  1  1;    % 7-9
              1 -1  0;   1  0  1;   0  1  1];   % 10-12
dirBeta   = [ 1  1  1;   1  1  1;   1  1  1;
             -1  1  1;  -1  1  1;  -1  1  1;
              1 -1  1;   1 -1  1;   1 -1  1;
              1  1 -1;   1  1 -1;   1  1 -1];
 
cAlpha = Miller(0,0,0,1,csAlpha,'hkil');    % alpha c-axis
aAlpha = Miller(1,1,-2,0,csAlpha,'UVTW');   % alpha a-axis [11-20]
 
% plane / direction of variant 1 (the reference)
p1 = Miller(planeBeta(1,1),planeBeta(1,2),planeBeta(1,3),csBeta,'hkl');
d1 = Miller(dirBeta(1,1),  dirBeta(1,2),  dirBeta(1,3),  csBeta,'uvw');
 
%% Reference OR = variant 1  (fit ONCE)
o1 = orientation.map(p1, cAlpha, d1, aAlpha);   % beta -> alpha, variant 1
 
%% Build all 12 variants IN TABLE ORDER
v = o1;  v = repmat(v,12,1);
for i = 1:12
    pI = Miller(planeBeta(i,1),planeBeta(i,2),planeBeta(i,3),csBeta,'hkl');
    dI = Miller(dirBeta(i,1),  dirBeta(i,2),  dirBeta(i,3),  csBeta,'uvw');
    % cubic operation carrying (p1,d1) -> (pI,dI): a genuine beta symmetry element
    s_i  = orientation.map(p1, pI, d1, dI);
    v(i) = o1 .* inv(s_i);            % variant i = reference composed with symmetry
end
 
%% Misorientation of each variant w.r.t. variant 1
ang = zeros(12,1); axStr = strings(12,1);
for i = 1:12
    m       = inv(v(1)) .* v(i);          % alpha-alpha disorientation
    ang(i)  = angle(m)./degree;
    a       = m.axis;                     % axis in variant-1 crystal frame
    axStr(i)= sprintf('[%.3f %.3f %.3f %.3f]', a.UVTW);
end
ang(1) = 0;  axStr(1) = "-";
 
T = table((1:12)', round(ang,2), axStr, ...
          'VariableNames', {'Variant','Angle_deg','Axis_UVTW'});
disp(T)
 
fprintf('\nDistinct angles (deg): ');
fprintf('%.2f  ', unique(round(ang(2:end),2)));  fprintf('\n');
fprintf('Expected: 10.53  60.00  60.83  63.26  90.00\n');

vv   = o1.variants;              % 12 variants, MTEX's own order
mm   = inv(vv(1)) .* vv;
angleangle = unique(round(angle(mm)./degree,2))';
 
% writetable(T,'burgers_misorientations.csv')   % uncomment to export