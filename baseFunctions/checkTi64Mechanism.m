function out = checkTi64Mechanism(eu, exp_para, options)
% checkTi64Mechanism  Verify the Ti64 alpha-lath DRP mechanism grain by grain.
%
%   out = checkTi64Mechanism(eu, exp_para)
%   out = checkTi64Mechanism(eu, exp_para, measured=..., phiOffset=..., ...)
%
% THE VERIFICATION TOOL for the model built by lathGeometryTi64.  For each
% KNOWN orientation it
%   (1) prints every reflecting feature with its geometry and amplitude, so
%       the physics can be checked number by number against section 2.3.2
%       (c-axis azimuth, the two stripe azimuths, the specular theta = 90-2Phi,
%       and which stripe should dominate);
%   (2) tiles the measured DRP next to the ISOLATED contribution of each
%       mechanism and the combined prediction, so the decomposition can be
%       checked by eye;
%   (3) optionally draws the predicted stripe centrelines on top of the
%       measured DRP.
%
% The columns are, in order:
%   measured | stripe PERP to c | stripe PARALLEL to c | plateau peaks | combined
% i.e. the "point 1" mechanism, the "point 4" mechanism, the R2-type specular
% points, and the full model.
%
% EBSD -> DRM REGISTRATION ------------------------------------------------
% The DRM azimuth zero and the EBSD sample X are generally mounted at some
% fixed rotation about the sample normal, and the chapter notes a residual
% rotational misalignment between the two measurements.  phiOffset adds a
% global azimuth offset to every prediction; sweep it until ONE value lines up
% ALL grains.  It does not affect the dictionary (indexing absorbs rotation
% about Z), only this comparison.
%
% Inputs
%   eu        K x 3 Bunge Euler angles [phi1 Phi phi2] in degrees, one per row.
%   exp_para  struct with the DRP grid, .crystal and .lath (lathGeometryTi64).
%
% Options
%   measured   K x 1 cell (or th_num x ph_num x K array) of measured DRPs.
%   phiOffset  global azimuth offset added to the predictions, deg (default 0).
%   overlay    true (default) draws the predicted stripe diameters on the
%              measured DRP.
%   scaleColumns  true (default) scales the three isolated columns by their
%              relative amplitude, so the tiles show which mechanism is
%              actually the bright one.  DRPsim_Ti64 max-normalises every
%              pattern it returns, so WITHOUT this the isolated columns would
%              all look equally bright and the intensity ratio - the whole
%              point of the mechanism - would be invisible.  Set false to get
%              each contribution at full contrast instead.
%   verbose    true (default) prints the per-grain feature table.
%   labels     K x 1 string row labels (default "grain k").
%
% Output
%   out  struct with fields
%          .reflectors  K x 1 cell of the lathReflectors struct arrays
%          .drp         K x 4 cell {perp, parallel, peaks, combined} DRPs,
%                       columns 1-3 scaled to a common brightness when
%                       scaleColumns is true
%          .stripeAz    K x 2 [azimuth of the PERP-to-c stripe, of the ||-to-c
%                       stripe] in degrees, phiOffset included
%          .stripeAmp   K x 2 matching amplitudes (before normalisation)
%          .bright      K x 4 [az th] of the predicted along-stripe MAXIMUM
%                       for the perp-to-c and the ||-to-c stripe.  This is the
%                       bevel asymmetry: the arc of a rounded crest is centred
%                       on the direction it faces, not on the DRP centre
%          .cAz         K x 1 c-axis azimuth in the DRP, phiOffset included
%
% See also lathGeometryTi64, lathReflectors, DRPsim_Ti64.
% -------------------------------------------------------------------------
arguments
    eu double
    exp_para struct
    options.measured = []
    options.phiOffset (1,1) double = 0
    options.overlay (1,1) logical = true
    options.scaleColumns (1,1) logical = true
    options.verbose (1,1) logical = true
    options.labels string = string([])
end

if size(eu,2) ~= 3
    error('checkTi64Mechanism:eu','eu must be K x 3 [phi1 Phi phi2] in degrees.');
end
if ~isfield(exp_para,'lath') || isempty(exp_para.lath)
    error('checkTi64Mechanism:notSetup', ...
        'exp_para.lath is missing. Run exp_para = lathGeometryTi64(exp_para) first.');
end
K = size(eu,1);

meas    = local_measuredToCell(options.measured, K);
hasMeas = ~isempty(meas);

if isempty(options.labels)
    labels = strings(K,1);
    for g = 1:K, labels(g) = sprintf('grain %d',g); end
else
    labels = options.labels(:);
end

% ---- four exp_para variants that isolate the mechanisms -----------------
% Each keeps exactly one group of features alive by zeroing the amplitudes of
% the others; the "combined" column is the untouched model.
epPerp = local_only(exp_para, "terraceRidge");   % point-1 mechanism
epPar  = local_only(exp_para, "basalEdge");      % point-4 mechanism
epPeak = local_only(exp_para, "plateaus");       % R2-type specular points
epAll  = exp_para;

colName = ["stripe PERP to c" "stripe || to c" "plateau peaks" "combined"];
epList  = {epPerp, epPar, epPeak, epAll};

out.reflectors = cell(K,1);
out.drp        = cell(K,4);
out.stripeAz   = nan(K,2);
out.stripeAmp  = zeros(K,2);
out.bright     = nan(K,4);
out.cAz        = nan(K,1);

for g = 1:K
    for jj = 1:4
        out.drp{g,jj} = double(DRPsim_Ti64(eu(g,1),eu(g,2),eu(g,3),epList{jj})) / 255;
    end
    R = lathReflectors(eu(g,1),eu(g,2),eu(g,3),exp_para.lath);
    out.reflectors{g} = R;

    % Put the three isolated columns back on a COMMON intensity scale.  Each
    % came out of DRPsim_Ti64 max-normalised, so without this the faint stripe
    % would look as bright as the dominant one.
    if options.scaleColumns
        grp = [local_groupAmp(R,"terraceRidge"), ...
               local_groupAmp(R,"basalEdge"), ...
               local_groupAmp(R,["basalPlateau","prismPlateau"])];
        ref = max(grp);
        if ref > 0
            for jj = 1:3
                out.drp{g,jj} = out.drp{g,jj} * (grp(jj)/ref);
            end
        end
    end

    % strongest feature of each stripe family (variant 1 only, for the table)
    [azT, ampT, brT] = local_strongest(R,"terraceRidge");
    [azB, ampB, brB] = local_strongest(R,"basalEdge");
    out.stripeAz(g,:)  = mod([azT azB] + options.phiOffset, 180);
    out.stripeAmp(g,:) = [ampT ampB];
    out.bright(g,:)    = [mod(brT(1) + options.phiOffset,360), brT(2), ...
                          mod(brB(1) + options.phiOffset,360), brB(2)];
    out.cAz(g)         = mod(local_cAz(R) + options.phiOffset, 360);

    if options.verbose
        local_printGrain(labels(g), eu(g,:), R, options.phiOffset);
    end
end

% ---- draw ---------------------------------------------------------------
nCol = hasMeas + 4;
fig  = figure('Name','checkTi64Mechanism', ...
    'Position',[60 60 min(1900,250*nCol) min(1050,250*K)]);
tl = tiledlayout(fig,K,nCol,'TileSpacing','compact','Padding','compact');

% DRP disk radii in the stereographic projection used by DRPdisp
Rout = cosd(exp_para.th_min) / (1 + sind(exp_para.th_min));
Rin  = cosd(exp_para.th_max) / (1 + sind(exp_para.th_max));

for g = 1:K
    col = 0;
    if hasMeas
        col = col + 1;
        nexttile(tl,(g-1)*nCol+col);
        DRPdisp(DRP_norm(meas{g}),exp_para);
        if options.overlay
            hold on
            % white  = stripe perpendicular to c ("point 1" mechanism)
            % yellow = stripe parallel to c      ("point 4" mechanism)
            % the circle marks where the BEVEL should put the stripe maximum
            local_diameter(out.stripeAz(g,1), Rin, Rout, 'w',  out.stripeAmp(g,1));
            local_diameter(out.stripeAz(g,2), Rin, Rout, 'y',  out.stripeAmp(g,2));
            local_brightMark(out.bright(g,1:2), exp_para, 'w', out.stripeAmp(g,1));
            local_brightMark(out.bright(g,3:4), exp_para, 'y', out.stripeAmp(g,2));
            hold off
        end
        if g == 1, title('measured','FontWeight','bold'); end
    end
    for jj = 1:4
        col = col + 1;
        nexttile(tl,(g-1)*nCol+col);
        DRPdisp(out.drp{g,jj},exp_para);
        if g == 1
            if jj == 4
                title(colName(jj),'FontWeight','bold');
            else
                title(colName(jj));
            end
        end
    end
end
title(tl, sprintf(['Ti64 lath mechanism  -  rows = orientations,  ', ...
    'phiOffset = %g deg   (overlay: white = perp. to c, yellow = || to c)'], ...
    options.phiOffset), 'FontWeight','bold');
end

% =========================================================================
function ep = local_only(ep, keep)
% Zero every amplitude except the requested group.
z = ["ampTerraceRidge" "ampBasalEdge" "ampBasalPlateau" "ampPrismPlateau"];
switch keep
    case "terraceRidge", alive = "ampTerraceRidge";
    case "basalEdge",    alive = "ampBasalEdge";
    case "plateaus",     alive = ["ampBasalPlateau" "ampPrismPlateau"];
    otherwise,           alive = z;
end
for f = z
    if ~any(strcmp(f, alive))
        ep.lath.(char(f)) = 0;
    end
end
% the pyramidal term is a peak contribution, so keep it only with the peaks
if ~any(strcmp("ampBasalPlateau", alive))
    ep.lath.pyramidFrac = 0;
end
end

% =========================================================================
function [az, amp, bright] = local_strongest(R, name)
% Azimuth, amplitude and bevel-brightest point (azimuth, theta) of the
% strongest feature of one family, variant 1.
az = NaN; amp = 0; bright = [NaN NaN];
for ii = 1:numel(R)
    if R(ii).name == name && R(ii).variant == 1 && R(ii).amp > amp
        amp    = R(ii).amp;
        az     = R(ii).stripeAz;
        bright = [R(ii).arcCentreAz, R(ii).arcCentreTh];
    end
end
end

% =========================================================================
function a = local_groupAmp(R, names)
% Brightest amplitude among a group of feature names (variant 1 only).
a = 0;
for ii = 1:numel(R)
    if any(R(ii).name == names) && R(ii).variant == 1
        a = max(a, R(ii).amp);
    end
end
end

% =========================================================================
function az = local_cAz(R)
% The c-axis azimuth is carried on the terraceRidge feature (its crest IS the
% in-surface c direction); fall back to NaN when that feature is absent.
az = NaN;
for ii = 1:numel(R)
    if R(ii).name == "terraceRidge" && R(ii).variant == 1
        az = R(ii).crestAz;
        return
    end
end
end

% =========================================================================
function local_diameter(az, Rin, Rout, colr, amp)
if isnan(az) || ~(amp > 0), return; end
lw = 0.8 + 2.2*min(amp,1);          % line thickness tracks predicted brightness
for s = [0 180]
    a = az + s;
    plot([Rin Rout]*cosd(a), [Rin Rout]*sind(a), '--', ...
        'Color', colr, 'LineWidth', lw);
end
end

% =========================================================================
function local_brightMark(bright, exp_para, colr, amp)
% Mark the predicted along-stripe maximum, in the same stereographic
% projection DRPdisp uses.  Skipped when the bevel centre falls outside the
% measured theta annulus (nothing to compare it against there).
az = bright(1); th = bright(2);
if isnan(az) || isnan(th) || ~(amp > 0), return; end
if th < exp_para.th_min || th > exp_para.th_max, return; end
r = cosd(th) / (1 + sind(th));
plot(r*cosd(az), r*sind(az), 'o', 'MarkerEdgeColor', colr, ...
    'LineWidth', 1.5, 'MarkerSize', 9);
end

% =========================================================================
function local_printGrain(label, eu, R, phiOffset)
fprintf('\n%s  [phi1 Phi phi2] = [%.1f %.1f %.1f]   (phiOffset = %g deg)\n', ...
    label, eu(1), eu(2), eu(3), phiOffset);
cAz = local_cAz(R);
if ~isnan(cAz)
    fprintf('  c-axis azimuth in the DRP: %.1f deg   (expect phi1 - 90 + offset)\n', ...
        mod(cAz + phiOffset,360));
end
% brightAz/brightTh are where the BEVEL puts the maximum along the stripe -
% the asymmetry the arc model introduces.  For a basalEdge, brightTh should
% come out at 90 - 2*Phi, the same place as the flat basal terrace's peak.
fprintf('  %-14s %-6s %9s %8s %8s %8s %8s %8s %8s  %s\n', ...
    'feature','kind','stripeAz','peakTh','expose','coinc','brightAz','brightTh','amp','');
for ii = 1:numel(R)
    f = R(ii);
    if isnan(f.arcCentreAz)
        brightAz = NaN;
    else
        brightAz = mod(f.arcCentreAz + phiOffset,360);
    end
    fprintf('  %-14s %-6s %9s %8s %8.3f %8s %8s %8s %8.4f  v%d\n', ...
        f.name, f.kind, ...
        local_num(mod(f.stripeAz + phiOffset,180)), local_num(f.peakTh), ...
        f.expose, local_num(f.coinc), ...
        local_num(brightAz), local_num(f.arcCentreTh), f.amp, f.variant);
end
[~,aT] = local_strongest(R,"terraceRidge");
[~,aB] = local_strongest(R,"basalEdge");
if aT > 0 || aB > 0
    if aT >= aB
        fprintf('  -> dominant stripe is PERPENDICULAR to c  (perp:par = %.2f : 1)\n', ...
            aT / max(aB,eps));
    else
        fprintf('  -> dominant stripe is PARALLEL to c       (par:perp = %.2f : 1)\n', ...
            aB / max(aT,eps));
    end
end
end

% =========================================================================
function s = local_num(x)
if isnan(x)
    s = '    -';
else
    s = sprintf('%8.1f',x);
end
end

% =========================================================================
function measured = local_measuredToCell(m, K)
if isempty(m)
    measured = {};
    return
end
if iscell(m)
    measured = m(:);
elseif ndims(m) == 3
    measured = squeeze(num2cell(m,[1 2]));
else
    measured = {m};
end
if numel(measured) ~= K
    error('checkTi64Mechanism:measured', ...
        'measured must supply one DRP per orientation (K = %d).',K);
end
end
