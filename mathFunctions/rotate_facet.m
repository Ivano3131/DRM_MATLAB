function [all_rot, faceW, pairW] = rotate_facet(Aa,Bb,Cc,facet)
% to generate the equivalent orientation vectors with respect to (Aa,Bb,Cc)
% euler angle rotation direction pair.
% input x is a faceting vector x = [u v w], in order u >= v >= w.
% output all_rot is the whole vector combination.
% input Euler angle triplets are in Bunge convention.
% Edit date: Aug 27, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
faceting_weights = [1 1 20]; %Ti64
faceting_weights = [1 1 1]; %Ti7

if nargin < 4
    error('faceting information is not determined')
end

u = facet(1); v = facet(2); w = facet(3); %weights take care of the shape

%Ti64--------------------------------------------------------------------

% remove the sides which are infinitely long?
vec_facets = [u v w; -u v w; u -v w; u v -w;...
    u w v; -u w v; u -w v; u w -v; v w u; -v w u; v -w u; v w -u;...
    v u w; -v u w; v -u w; v u -w; w u v; -w u v; w -u v; w u -v;...
    w v u; -w v u; w -v u; w v -u];
vec_facets_shortened = [u v w; -u v w; u -v w; u v -w;...
    u w v; -u w v; u -w v; u w -v; v w u; -v w u; v -w u; v w -u;...
    w v u; -w v u; w -v u; w v -u];
vec_facets = vec_facets_shortened;
vec_eq = unique(vec_facets,'rows');
%vec_eq
[~, ia, ~] = unique(vec_facets, 'rows');
%ia
% this is for cubic system
face_u_dir = faceting_weights(2)*faceting_weights(3);
face_v_dir = faceting_weights(1)*faceting_weights(3);
face_w_dir = faceting_weights(1)*faceting_weights(2);
faceW = [face_u_dir; face_u_dir; face_u_dir; face_u_dir;...
    face_u_dir; face_u_dir; face_u_dir; face_u_dir; face_w_dir; face_w_dir; face_w_dir; face_w_dir;
    face_v_dir; face_v_dir; face_v_dir; face_v_dir; face_v_dir; face_v_dir; face_v_dir; face_v_dir;...
   face_w_dir; face_w_dir; face_w_dir; face_w_dir];
faceW_shortened = [face_u_dir; face_u_dir; face_u_dir; face_u_dir;...
    face_u_dir; face_u_dir; face_u_dir; face_u_dir; face_w_dir; face_w_dir; face_w_dir; face_w_dir;
   face_w_dir; face_w_dir; face_w_dir; face_w_dir];
faceW = faceW_shortened;
faceW = faceW(ia,:);
pairW = [faceting_weights(1);faceting_weights(1);faceting_weights(1);faceting_weights(1);...
    faceting_weights(1);faceting_weights(1);faceting_weights(1);faceting_weights(1); faceting_weights(3); faceting_weights(3);...
    faceting_weights(3); faceting_weights(3); faceting_weights(2);faceting_weights(2);faceting_weights(2);faceting_weights(2);faceting_weights(2);faceting_weights(2);faceting_weights(2);faceting_weights(2);
    faceting_weights(3);faceting_weights(3);faceting_weights(3);faceting_weights(3)];
pairW_shortened = [faceting_weights(1);faceting_weights(1);faceting_weights(1);faceting_weights(1);...
    faceting_weights(1);faceting_weights(1);faceting_weights(1);faceting_weights(1); faceting_weights(3); faceting_weights(3);...
    faceting_weights(3); faceting_weights(3); faceting_weights(3);faceting_weights(3);faceting_weights(3);faceting_weights(3)];
pairW = pairW_shortened;
pairW = pairW(ia,:);

%%cuboid with weights
%u = facet(1); v = facet(2); w = facet(3);
%vec_eq = unique([u v w; u v -w; -u -v w], 'rows');

%faceW = [faceting_weights(2)*faceting_weights(3), faceting_weights(2)*faceting_weights(3), faceting_weights(1)*faceting_weights(3), ...
%    faceting_weights(1)*faceting_weights(3),faceting_weights(1)*faceting_weights(2),faceting_weights(1)*faceting_weights(2)];

%pairW = zeros(3, 3); % Initialize pairW matrix

%for i = 1:3
%    ai = axis_id(vec_eq(i, :));
%    for j = i+1:3
%        aj = axis_id(vec_eq(j,:));
%        if ai == aj
%            continue;
%        end
%        missingAxis = 6 - ai - aj;
%        if missingAxis == 1
%            w = faceting_weights(1);
%        elseif missingAxis == 2
%            w = faceting_weights(2);
%        else
%            w = faceting_weights(3);
%        end
%        pairW(i,j) = w;
%        pairW(j,i) = w;
%    end
%end

%Ti64 ---------------------------------------------------------------

%Ti7
vec_facets = [w w u; w w -u; u -u w; -u u w; 2*u u w; -2*u -u w; u 2*u w; -u -2*u w];
vec_eq = unique(vec_facets, 'rows');
faceW = ones(size(vec_eq,1),1);
pairW = ones(size(vec_eq,1),1);

faceW = faceW / sum(faceW);
if any(pairW(:) > 0)
    pairW = pairW / sum(pairW(:));
end

N_vec = length(vec_eq);
%% all_rot=EulerRotate(all11x,Aa,Bb,Cc);
all_rot = normr(EulerRotate(vec_eq,Aa,Bb,Cc));
for ii = 1:N_vec
    if all_rot(ii,3)<0
        all_rot(ii,:)=-all_rot(ii,:);
    end
end

end

%{
function ax = axis_id(v)
[~, ax] = max(abs(v));
end
%}

%u_hcp = facet(1); v_hcp = facet(2); w_hcp = facet(3); r_hcp = facet(4);

%base = facet; %treats it as if it was a cube
%P = base(perms(1:4), :);
%S = dec2bin(0:2^4-1)-'0';
%S = 1 - 2*S;

%vec_eq_hcp = unique(reshape(permute(P, [1 3 2]) .* permute(S, [3 1 2]), [], 4), 'rows');


%% this is probably THE place to create the terrace features

%R = [
%    u, v, t, w;
%    -t, u, v, w;
%    -v, -t, u, w;
%    -u, -v, -t, w;
%    t, -u, -v, w;
%    v, t, -u, w;
%    ];

%M = R(:, [1 3 2 4]); %basal plane
%ops = [R; M]; %antipodal
%ops = [ops; -ops];
%vec_eq_hcp = unique(ops, 'rows');

% plate-like terrace features
%R = [u, v, t, w;
%    -u, -v, -t, w;
%    u, v, t, -w];
%vec_eq_terrace = unique(R, 'rows');
%vec_eq_hcp = vec_eq_terrace;

%N_vec = length(vec_eq_hcp);
%all_rot = normr(EulerRotate(vec_eq_hcp,Aa, BB, Cc));
%for ii = 1:N_vec
%    if all_rot(ii, 3) < 0
%        all_rot(ii, :) = -all_rot(ii, :);
%    end

%end