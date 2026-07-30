function out = checkCubeMechanism(eu, exp_para, options)
% checkCubeMechanism  Verify the Ti64 cube-lath mechanism grain by grain.
%
%   out = checkCubeMechanism(eu, exp_para)
%   out = checkCubeMechanism(eu, exp_para, measured=..., phiOffset=..., ...)
%
% THE VERIFICATION TOOL for the model built by cubeGeometryTi64.  For each
% KNOWN orientation it
%   (1) prints every edge with its exposure and geometry, so the etch cut can
%       be checked number by number (c-axis azimuth, edge tilt, stripe
%       azimuth, exposed length fraction, and which family should dominate);
%   (2) tiles the measured DRP next to the ISOLATED contribution of each edge
%       family and the combined prediction;
%   (3) optionally draws the predicted stripe centrelines on the measured DRP.
%
% The columns are, in order:
%   measured | stripe PERP to c | stripe PARALLEL to c | all edges
%
% where the second column is the four edges PARALLEL to c (broad face ^ end
% face) - the Fig 2.12 point-1 mechanism - and the third is the eight edges
% lying in the basal plane (the basal rim) - the point-4 mechanism.
%
% WHAT TO EXPECT ----------------------------------------------------------
% Points 1-3 have Phi >= 72 deg: c lies near the surface, the box stands on
% edge, so the edges parallel to c sit at the top of the etch cut and the
% PERPENDICULAR-to-c stripe should dominate.  Point 4 has Phi = 37.8 deg: the
% basal face is near-horizontal and forms the top, so its rim is exposed and
% the PARALLEL-to-c stripe should dominate.  Nothing in the model asserts
% this - it comes out of the printed `expose` column.
%
% EBSD -> DRM REGISTRATION ------------------------------------------------
% The DRM azimuth zero and the EBSD sample X are generally mounted at some
% fixed rotation about the sample normal.  phiOffset adds a global azimuth
% offset to every prediction; sweep it until ONE value lines up ALL grains.
% It does not affect indexing (matchDRPcube absorbs rotation about Z), only
% this comparison.
%
% Inputs
%   eu        K x 3 Bunge Euler angles [phi1 Phi phi2] in degrees, one per row.
%   exp_para  struct with the DRP grid, .crystal and .cube (cubeGeometryTi64).
%
% Options
%   measured   K x 1 cell (or th_num x ph_num x K array) of measured DRPs.
%   phiOffset  global azimuth offset added to the predictions, deg (default 0).
%   overlay    true (default) draws the predicted stripe diameters on the
%              measured DRP.
%   scaleColumns  true (default) scales the two isolated columns by their
%              relative amplitude, so the tiles show which family is actually
%              the bright one.  DRPsim_cube max-normalises every pattern it
%              returns, so WITHOUT this both columns would look equally
%              bright and the intensity ratio - the whole point - would be
%              invisible.
%   verbose    true (default) prints the per-grain edge table.
%   labels     K x 1 string row labels (default "grain k").
%
% Output
%   out  struct with fields
%          .edges      K x 1 cell of the cubeReflectors struct arrays
%          .drp        K x 3 cell {perp-to-c, parallel-to-c, all}, columns 1-2
%                      scaled to a common brightness when scaleColumns is true
%          .stripeAz   K x 2 [azimuth of the PERP-to-c stripe, of the ||-to-c
%                      stripe] in degrees, phiOffset included
%          .stripeAmp  K x 2 matching amplitudes (before normalisation)
%          .expose     K x 2 exposed length fraction of the two families
%          .bright     K x 4 [az th az th] of the predicted along-stripe
%                      MAXIMUM for each family - the bevel asymmetry
%          .cAz        K x 1 c-axis azimuth in the DRP, phiOffset included
%
% See also cubeGeometryTi64, cubeReflectors, DRPsim_cube, matchDRPcube.
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
    error('checkCubeMechanism:eu','eu must be K x 3 [phi1 Phi phi2] in degrees.');
end
if ~isfield(exp_para,'cube') || isempty(exp_para.cube)
    error('checkCubeMechanism:notSetup', ...
        'exp_para.cube is missing. Run exp_para = cubeGeometryTi64(exp_para) first.');
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

% ---- three exp_para variants that isolate the edge families -------------
% Isolation is by AXIS, not by amplitude knob: axis 1 is parallel to c, axes
% 2 and 3 lie in the basal plane.
epPerp = local_onlyAxes(exp_para, 1);       % edges || c   -> stripe PERP to c
epPar  = local_onlyAxes(exp_para, [2 3]);   % basal rim    -> stripe || to c
epAll  = exp_para;

colName = ["stripe PERP to c" "stripe || to c" "all edges"];
epList  = {epPerp, epPar, epAll};

out.edges     = cell(K,1);
out.drp       = cell(K,3);
out.stripeAz  = nan(K,2);
out.stripeAmp = zeros(K,2);
out.expose    = zeros(K,2);
out.bright    = nan(K,4);
out.cAz       = nan(K,1);

for g = 1:K
    for jj = 1:3
        out.drp{g,jj} = double(DRPsim_cube(eu(g,1),eu(g,2),eu(g,3),epList{jj})) / 255;
    end
    E = cubeReflectors(eu(g,1),eu(g,2),eu(g,3),exp_para.cube);
    out.edges{g} = E;

    [azT, ampT, brT, expT] = local_strongest(E, 1);
    [azB, ampB, brB, expB] = local_strongest(E, [2 3]);

    % Put the two isolated columns back on a COMMON intensity scale.
    if options.scaleColumns
        ref = max([ampT ampB]);
        if ref > 0
            out.drp{g,1} = out.drp{g,1} * (ampT/ref);
            out.drp{g,2} = out.drp{g,2} * (ampB/ref);
        end
    end

    out.stripeAz(g,:)  = mod([azT azB] + options.phiOffset, 180);
    out.stripeAmp(g,:) = [ampT ampB];
    out.expose(g,:)    = [expT expB];
    out.bright(g,:)    = [mod(brT(1) + options.phiOffset,360), brT(2), ...
                          mod(brB(1) + options.phiOffset,360), brB(2)];
    out.cAz(g)         = mod(E(1).cAz + options.phiOffset, 360);

    if options.verbose
        local_printGrain(labels(g), eu(g,:), E, options.phiOffset);
    end
end

% ---- draw ---------------------------------------------------------------
nCol = hasMeas + 3;
fig  = figure('Name','checkCubeMechanism', ...
    'Position',[60 60 min(1900,260*nCol) min(1050,260*K)]);
tl = tiledlayout(fig,K,nCol,'TileSpacing','compact','Padding','compact');

% DRP disk radii in the stereographic projection used by DRPdisp
Rout = cosd(exp_para.th_min) / (1 + sind(exp_para.th_min));   % outer = th_min
Rin  = cosd(exp_para.th_max) / (1 + sind(exp_para.th_max));   % inner = th_max

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
            % the circle marks where the bevel should put the stripe maximum
            local_diameter(out.stripeAz(g,1), Rin, Rout, 'w', out.stripeAmp(g,1), exp_para);
            local_diameter(out.stripeAz(g,2), Rin, Rout, 'y', out.stripeAmp(g,2), exp_para);
            local_brightMark(out.bright(g,1:2), exp_para, 'w', out.stripeAmp(g,1));
            local_brightMark(out.bright(g,3:4), exp_para, 'y', out.stripeAmp(g,2));
            hold off
        end
        if g == 1, title('measured','FontWeight','bold'); end
    end
    for jj = 1:3
        col = col + 1;
        nexttile(tl,(g-1)*nCol+col);
        DRPdisp(out.drp{g,jj},exp_para);
        if g == 1
            if jj == 3
                title(colName(jj),'FontWeight','bold');
            else
                title(colName(jj));
            end
        end
    end
end
title(tl, sprintf(['Ti64 cube-lath mechanism  -  rows = orientations,  ', ...
    'etchFrac = %g,  phiOffset = %g deg   ', ...
    '(overlay: white = perp. to c, yellow = || to c)'], ...
    exp_para.cube.etchFrac, options.phiOffset), 'FontWeight','bold');
end

% =========================================================================
function ep = local_onlyAxes(ep, keepAxes)
% Silence every edge whose axis is not in keepAxes.
kill = ~ismember(ep.cube.edgeAxis, keepAxes);
ep.cube.edgeAmp(kill) = 0;
end

% =========================================================================
function [az, amp, bright, expose] = local_strongest(E, axesWanted)
% Azimuth, amplitude, bevel-brightest point and exposure of the strongest
% edge among the requested axes.
az = NaN; amp = 0; bright = [NaN NaN]; expose = 0;
for ii = 1:numel(E)
    if ismember(E(ii).axisId, axesWanted) && E(ii).amp > amp
        amp    = E(ii).amp;
        az     = E(ii).stripeAz;
        bright = [E(ii).arcCentreAz, E(ii).arcCentreTh];
        expose = E(ii).expose;
    end
end
end

% =========================================================================
function a = local_dispAz(az, exp_para)
% DRPdisp.m:29 lays its ph_num columns out over 0:360/(ph_num-1):360, while
% the data grid is 0:360/ph_num:360-360/ph_num.  Overlays therefore have to
% be stretched by the same factor or they drift ~3 deg from the pattern they
% are meant to mark by the time the azimuth wraps.
a = az * exp_para.ph_num / (exp_para.ph_num - 1);
end

% =========================================================================
function local_diameter(az, Rin, Rout, colr, amp, exp_para)
if isnan(az) || ~(amp > 0), return; end
lw = 0.8 + 2.2*min(amp,1);          % line thickness tracks predicted brightness
for s = [0 180]
    a = local_dispAz(az + s, exp_para);
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
a = local_dispAz(az, exp_para);
plot(r*cosd(a), r*sind(a), 'o', 'MarkerEdgeColor', colr, ...
    'LineWidth', 1.5, 'MarkerSize', 9);
end

% =========================================================================
function local_printGrain(label, eu, E, phiOffset)
fprintf('\n%s  [phi1 Phi phi2] = [%.1f %.1f %.1f]   (phiOffset = %g deg)\n', ...
    label, eu(1), eu(2), eu(3), phiOffset);
fprintf('  c-axis azimuth in the DRP: %.1f deg   (expect phi1 - 90 + offset)\n', ...
    mod(E(1).cAz + phiOffset, 360));
fprintf('  c-axis polar angle Phi:    %.1f deg\n', E(1).Phi);
fprintf('  %-12s %5s %7s %9s %8s %9s %9s %8s\n', ...
    'edge','axis','tilt','stripeAz','expose','brightAz','brightTh','amp');
for ii = 1:numel(E)
    f = E(ii);
    fprintf('  %-12s %5d %7.1f %9.1f %8.3f %9s %9s %8.4f\n', ...
        f.family, f.axisId, f.tilt, mod(f.stripeAz + phiOffset,180), ...
        f.expose, local_num(mod(f.arcCentreAz + phiOffset,360)), ...
        local_num(f.arcCentreTh), f.amp);
end
[~,aT] = local_strongest(E, 1);
[~,aB] = local_strongest(E, [2 3]);
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
    s = '        -';
else
    s = sprintf('%9.1f',x);
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
    error('checkCubeMechanism:measured', ...
        'measured must supply one DRP per orientation (K = %d).',K);
end
end
