function simDRP = DRPsim_hcp(eu1,eu2,eu3,exp_para)
% DRPsim_hcp  Simulate a Directional Reflectance Pattern for an HCP grain.
%
%   simDRP = DRPsim_hcp(eu1,eu2,eu3,exp_para)
%
% Two selectable surface-topology models (exp_para.reflectorModel):
%
%   "facet" (default) - flat crystallographic facets.  Builds (a) specular
%       "peaks" where a facet mirrors light into the fixed overhead camera, and
%       (b) great-circle "bands" from the edges between pairs of facets.  Facet
%       normals come from rotate_facet_hcp / facetNormalsHCP.
%
%   "ridge" - smooth ELONGATED RIDGES (half-pipe / cylindrical topology, as seen
%       by AFM/SEM).  A cylinder with axis t has normals covering the whole great
%       circle perpendicular to t, so it reflects the camera along that great
%       circle -> a STRIPE in the DRP whose pole is the ridge axis.  Axes come
%       from ridgeAxesHCP (exp_para.ridgeUVTW).  When the ridge lies in the
%       surface, the stripe projects to a near-diameter line through the DRP
%       centre.
%
%   "both" - facet peaks/bands AND ridge stripes, combined by max.
%
% DRP geometry: each DRP cell (theta,phi) maps to the probe direction
% thph2vec(45 + theta/2, phi), i.e. the microfacet normal that would send the
% illumination specularly into the camera.
%
% Relevant exp_para fields
%   th_max, th_min, th_num, ph_num   DRP grid (illumination sampling)
%   fitting_para = [i_Main, i_facet, sd_Main, sd_facet, ...]
%       i_Main   peak amplitude,  sd_Main  peak half-width (deg, "specular")
%       i_facet  band amplitude,  sd_facet band half-width (deg)
%   facetNormals, faceW, pairW       from facetNormalsHCP   (facet model)
%   ridgePlaneNormals, ridgeW        from ridgePlanesHCP    (ridge planeTrace)
%   ridgeAxes, ridgeW                from ridgeAxesHCP      (ridge direction)
%   peakModel  "specular" (default, Cauchy lobe) or "cosine" (broad cos lobe)
%   reflectorModel "facet" (default) | "ridge" | "both"
%   ridgeMode  "planeTrace" (default) ridge = surface trace m x Z of a crystal
%                plane (stripe ALWAYS through DRP centre); or "direction" ridge
%                = a rigidly-rotated crystallographic direction (off-centre).
%   ridgeFitting = [i_stripe, sd_stripe]  stripe amplitude / cross-width (deg);
%                                         default [1, 6]
%   ridgeTaper   >=0 exponent that fades the stripe toward the DRP edge (low
%                theta); 0 = uniform stripe (default), larger = more central
% -------------------------------------------------------------------------
arguments
    eu1 double
    eu2 double
    eu3 double
    exp_para struct
end

th_max = exp_para.th_max;
th_min = exp_para.th_min;
th_num = exp_para.th_num;
ph_num = exp_para.ph_num;

i_Main   = exp_para.fitting_para(1);
i_facet  = exp_para.fitting_para(2);
sd_Main  = exp_para.fitting_para(3);
sd_facet = exp_para.fitting_para(4);

if isfield(exp_para,'peakModel') && ~isempty(exp_para.peakModel)
    peakModel = string(exp_para.peakModel);
else
    peakModel = "specular";
end
if isfield(exp_para,'reflectorModel') && ~isempty(exp_para.reflectorModel)
    reflectorModel = string(exp_para.reflectorModel);
else
    reflectorModel = "facet";
end
includeFacet = reflectorModel == "facet" || reflectorModel == "both";
includeRidge = reflectorModel == "ridge" || reflectorModel == "both";

% facet normals in the specimen frame for this orientation (facet model only)
if includeFacet
    [rot_facet, faceW, pairW] = rotate_facet_hcp(eu1,eu2,eu3,exp_para);
end

th_step = (th_max - th_min) / (th_num - 1);
ph_step = 360 / ph_num;

cauchy = @(p,x) p(1) ./ (1 + (x./p(2)).^2);   % i / (1 + (dist/sd)^2)

% DRP sampling grid (one extra leading theta row, cropped off at the end)
th_range = th_min-th_step : th_step : ...
    (floor((90-th_min)/th_step)-1)*th_step + th_min;
th_DRP = repmat(transpose(th_range),1,ph_num);
ph_DRP = repmat(0:ph_step:360-ph_step,length(th_range),1);
simDRP = zeros(length(th_range),ph_num);
vec_DRP = zeros(3,length(th_range),ph_num);

% probe direction for every DRP cell = facet normal that specularly reflects
for ii = 1:length(th_range)
    for jj = 1:ph_num
        tmp_vec = thph2vec(45 + th_DRP(ii,jj)/2, ph_DRP(ii,jj));
        vec_DRP(:,ii,jj) = normr(tmp_vec);
    end
end
vx = squeeze(vec_DRP(1,:,:));
vy = squeeze(vec_DRP(2,:,:));
vz = squeeze(vec_DRP(3,:,:));

% ==== FACET model: specular peaks + facet-pair bands ====================
if includeFacet
    % ---- major reflectance peaks --------------------------------------
    for ii = 1:size(rot_facet,1)
        n = rot_facet(ii,:);
        cosTheta = n(1)*vx + n(2)*vy + n(3)*vz;     % alignment with probe
        cosTheta = max(min(cosTheta,1),-1);
        switch peakModel
            case "cosine"
                temp = faceW(ii) * i_Main * max(cosTheta,0);
            otherwise   % "specular": sharp Cauchy lobe of half-width sd_Main
                peakDist = acosd(cosTheta);         % deg from facet normal
                temp = faceW(ii) * cauchy([i_Main, sd_Main], peakDist);
        end
        simDRP = max(simDRP, temp);
    end

    % ---- great-circle bands between facet pairs -----------------------
    for ii = 1:size(rot_facet,1)
        for jj = ii+1:size(rot_facet,1)
            vec_1 = rot_facet(ii,:);
            vec_2 = rot_facet(jj,:);
            if norm(vec_1-vec_2) < 1e-3 || norm(vec_1+vec_2) < 1e-3
                continue
            end
            gcnorm = normr(cross(vec_1,vec_2));
            bandCos = max(min(gcnorm(1)*vx + gcnorm(2)*vy + gcnorm(3)*vz,1),-1);
            bandDist = asind(bandCos);              % deg from great circle
            w = max(pairW(ii),pairW(jj));
            temp = w * cauchy([i_facet, sd_facet], bandDist);
            simDRP = max(simDRP, temp);
        end
    end
end

% ==== RIDGE model: cylindrical ridges -> great-circle stripes ===========
if includeRidge
    if isfield(exp_para,'ridgeFitting') && ~isempty(exp_para.ridgeFitting)
        i_stripe  = exp_para.ridgeFitting(1);
        sd_stripe = exp_para.ridgeFitting(2);
    else
        i_stripe = 1; sd_stripe = 6;
    end
    if isfield(exp_para,'ridgeTaper') && ~isempty(exp_para.ridgeTaper)
        ridgeTaper = exp_para.ridgeTaper;
    else
        ridgeTaper = 0;
    end
    if isfield(exp_para,'ridgeMode') && ~isempty(exp_para.ridgeMode)
        ridgeMode = string(exp_para.ridgeMode);
    else
        ridgeMode = "planeTrace";
    end

    % assemble the specimen-frame ridge axes tList (each row a unit line)
    if ridgeMode == "direction"
        % axis IS a crystallographic direction (rigidly rotated); the stripe
        % need not pass through the DRP centre.
        if ~isfield(exp_para,'ridgeAxes') || isempty(exp_para.ridgeAxes)
            error('DRPsim_hcp:noRidgeAxes', ...
                'ridgeMode="direction" needs exp_para.ridgeAxes (run ridgeAxesHCP).');
        end
        tList = normr(EulerRotate(exp_para.ridgeAxes, eu1, eu2, eu3));
    else
        % planeTrace (default): ridge axis = surface trace of a crystal plane,
        % t = m x Z.  This is ALWAYS horizontal -> stripe through the centre.
        if ~isfield(exp_para,'ridgePlaneNormals') || isempty(exp_para.ridgePlaneNormals)
            error('DRPsim_hcp:noRidgePlanes', ...
                'ridgeMode="planeTrace" needs exp_para.ridgePlaneNormals (run ridgePlanesHCP).');
        end
        mS = normr(EulerRotate(exp_para.ridgePlaneNormals, eu1, eu2, eu3));
        tList = cross(mS, repmat([0 0 1],size(mS,1),1), 2);   % m x Z, in-plane
    end
    ridgeW = exp_para.ridgeW;

    taper = 1;
    if ridgeTaper > 0
        taper = max(vz,0).^ridgeTaper;   % fade the stripe toward the DRP edge
    end
    for ii = 1:size(tList,1)
        t = tList(ii,:);
        nt = norm(t);
        if nt < 1e-6            % plane parallel to surface -> no ridge line
            continue
        end
        t = t / nt;
        stripeCos = max(min(t(1)*vx + t(2)*vy + t(3)*vz,1),-1);
        stripeDist = asind(stripeCos);              % deg from the great circle
        temp = ridgeW(ii) * cauchy([i_stripe, sd_stripe], stripeDist) .* taper;
        simDRP = max(simDRP, temp);
    end
end

% normalise to [0,1], crop the leading theta row, return as uint8 image
simDRP = simDRP - min(simDRP(:));
if max(simDRP(:)) > 0
    simDRP = simDRP / max(simDRP(:));
end
simDRP = uint8(simDRP(2:th_num+1,:) * 255);
end
