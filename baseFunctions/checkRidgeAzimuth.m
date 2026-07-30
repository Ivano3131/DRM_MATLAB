function info = checkRidgeAzimuth(eu, measured, exp_para, options)
% checkRidgeAzimuth  Overlay the predicted ridge-stripe centreline(s) on the
%                    measured DRPs, to confirm the ridge plane and the
%                    EBSD -> DRM azimuth registration.
%
%   info = checkRidgeAzimuth(eu, measured, exp_para)
%   info = checkRidgeAzimuth(eu, measured, exp_para, ridgeHKIL=..., phiOffset=...)
%
% For the plane-trace ridge model the stripe is a diameter through the DRP
% centre at azimuth phi_s where the surface trace t = m x Z is perpendicular to
% the probe.  This draws that predicted diameter on top of each measured DRP so
% you can see whether it lies along the real streak.
%
% Because the DRM azimuth zero and the EBSD sample-X are generally mounted at a
% fixed rotation about the sample normal, sweep `phiOffset` until ONE value
% lines up ALL grains at once: that value is your EBSD -> DRM registration, and
% a consistent fit confirms both the ridge plane and the orientation convention.
% (A single plane like basal gives one line; prism/pyramidal give several.)
%
% Inputs
%   eu        K x 3 Bunge Euler angles [phi1 Phi phi2] in degrees (from EBSD).
%   measured  K x 1 cell (or th_num x ph_num x K array) of measured DRPs.
%   exp_para  struct (uses .crystal, .th_min, .th_max, and .ridgeHKIL if set).
%
% Options
%   ridgeHKIL  1 x 4 (or M x 4) {h k i l} ridge plane (default: exp_para.ridgeHKIL,
%              else basal [0 0 0 1]).
%   phiOffset  global azimuth offset added to every predicted line (deg).
%   lineColor  overlay colour (default 'w').
%
% Output
%   info  K x 1 struct with .eu and .predAz (predicted stripe azimuths, deg).
% -------------------------------------------------------------------------
arguments
    eu double
    measured
    exp_para struct
    options.ridgeHKIL double = []
    options.phiOffset (1,1) double = 0
    options.lineColor = 'w'
end

K = size(eu,1);
meas = local_measuredToCell(measured, K);

if isempty(options.ridgeHKIL)
    if isfield(exp_para,'ridgeHKIL') && ~isempty(exp_para.ridgeHKIL)
        hkil = exp_para.ridgeHKIL;
    else
        hkil = [0 0 0 1];
    end
else
    hkil = options.ridgeHKIL;
end

ep = exp_para; ep.ridgeHKIL = hkil; ep = ridgePlanesHCP(ep);
Nrm = ep.ridgePlaneNormals;                 % crystal-frame plane normals

% DRP disk radii (stereographic, as in DRPdisp): outer = th_min, inner = th_max
Rout = cosd(exp_para.th_min) / (1 + sind(exp_para.th_min));
Rin  = cosd(exp_para.th_max) / (1 + sind(exp_para.th_max));

fig = figure('Name','checkRidgeAzimuth', ...
    'Position',[100 100 min(1800,340*K) 380]); %#ok<NASGU>
tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');

info = struct('eu',cell(K,1),'predAz',cell(K,1));
for g = 1:K
    nexttile(tl);
    DRPdisp(DRP_norm(meas{g}), exp_para);
    hold on;

    mS = normr(EulerRotate(Nrm, eu(g,1), eu(g,2), eu(g,3)));  % specimen frame
    azs = [];
    for ii = 1:size(mS,1)
        t = cross(mS(ii,:), [0 0 1]);        % surface trace, m x Z
        nt = norm(t);
        if nt < 1e-6, continue; end          % plane parallel to surface
        t = t / nt;
        phi_s = mod(atan2d(-t(1), t(2)) + options.phiOffset, 360);
        azs(end+1) = mod(phi_s,180); %#ok<AGROW>
        for s = [0 180]
            a = phi_s + s;
            plot([Rin Rout]*cosd(a), [Rin Rout]*sind(a), '--', ...
                'Color', options.lineColor, 'LineWidth', 1.5);
        end
    end
    hold off;
    title(sprintf('grain %d: pred az = %s^o', g, mat2str(round(unique(azs),1))));
    info(g).eu = eu(g,:);
    info(g).predAz = azs;
end
sgtitle(sprintf('predicted {%d %d %d %d} ridge stripe on measured DRP  (phiOffset = %g^o)', ...
    hkil(1,1),hkil(1,2),hkil(1,3),hkil(1,4),options.phiOffset));
end

% -------------------------------------------------------------------------
function measured = local_measuredToCell(m, K)
if iscell(m)
    measured = m(:);
elseif ndims(m) == 3
    measured = squeeze(num2cell(m,[1 2]));
else
    measured = {m};
end
if numel(measured) ~= K
    error('checkRidgeAzimuth:measured', ...
        'measured must supply one DRP per orientation (K = %d).',K);
end
end
