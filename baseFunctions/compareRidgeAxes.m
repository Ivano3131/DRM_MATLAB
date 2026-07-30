function results = compareRidgeAxes(eu, exp_para, options)
% compareRidgeAxes  A/B simulated DRPs across candidate ridge families.
%
%   results = compareRidgeAxes(eu, exp_para)
%   results = compareRidgeAxes(eu, exp_para, families=..., labels=..., ...
%                              measured=..., mode=..., ridgeFitting=..., ridgeTaper=...)
%
% Ridge-topology counterpart of compareFacetFamilies.  For one or more KNOWN
% orientations it simulates the DRP with the "ridge" reflector model and tiles
% the candidates next to the measured DRP, so the ridge crystallography that
% reproduces the observed streak can be picked by eye.
%
% Two ridge geometries (options.mode):
%   "planeTrace" (default) - ridges are the SURFACE TRACES of a crystallographic
%       plane, axis t = m x Z.  t is horizontal for every orientation, so the
%       stripe ALWAYS passes through the DRP centre (matches real surface
%       ridges).  families are {h k i l} PLANES.  Default candidates: basal,
%       prism {10-10}, pyramidal {10-11}, {10-12}.
%   "direction" - ridge axis is a rigidly-rotated crystallographic DIRECTION
%       <u v t w>; the stripe generally sits off-centre.  families are <uvtw>.
%
% Inputs
%   eu        K x 3 Bunge Euler angles [phi1 Phi phi2] in degrees (one per row).
%   exp_para  struct with the DRP grid + fitting fields and (ideally) .crystal.
%
% Options
%   mode          "planeTrace" (default) | "direction"
%   families      1 x N cell of candidate families (M x 4 each); default depends
%                 on mode (planes for planeTrace, directions for direction).
%   labels        1 x N string column titles (default derived from families).
%   measured      K x 1 cell (or th x ph x K array) of measured DRPs (optional).
%   ridgeFitting  [i_stripe, sd_stripe] stripe amplitude / cross-width (deg).
%   ridgeTaper    >=0 fade of the stripe toward the DRP edge (default 0).
%
% Output
%   results   struct array (K x N) with fields .eu .label .family .drp .nAxes
% -------------------------------------------------------------------------
arguments
    eu double
    exp_para struct
    options.mode string = "planeTrace"
    options.families cell = {}
    options.labels string = string([])
    options.measured = []
    options.ridgeFitting double = []
    options.ridgeTaper double = []
end

if size(eu,2) ~= 3
    error('compareRidgeAxes:eu','eu must be K x 3 [phi1 Phi phi2] in degrees.');
end
K = size(eu,1);
mode = options.mode;

% default candidate list depends on the geometry
families = options.families;
if isempty(families)
    if mode == "direction"
        families = { [1 1 -2 0], [1 0 -1 0], [0 0 0 1], [1 1 -2 3] };
    else
        families = { [0 0 0 1], [1 0 -1 0], [1 0 -1 1], [1 0 -1 2] };
    end
end
N = numel(families);

if isempty(options.labels)
    labels = strings(1,N);
    for jj = 1:N
        labels(jj) = local_familyLabel(families{jj}, mode);
    end
else
    labels = options.labels;
end

measured = local_measuredToCell(options.measured, K);
hasMeas = ~isempty(measured);

% force the ridge model / geometry for this comparison
exp_para.reflectorModel = "ridge";
exp_para.ridgeMode = mode;
if ~isempty(options.ridgeFitting), exp_para.ridgeFitting = options.ridgeFitting; end
if ~isempty(options.ridgeTaper),   exp_para.ridgeTaper   = options.ridgeTaper;   end

results = struct('eu',[],'label',[],'family',[],'drp',[],'nAxes',[]);
for jj = 1:N
    ep = exp_para;
    if mode == "direction"
        ep.ridgeUVTW = families{jj};  ep.ridgeWeights = [];
        ep = ridgeAxesHCP(ep);
        nAx = size(ep.ridgeAxes,1);
    else
        ep.ridgeHKIL = families{jj};  ep.ridgeWeights = [];
        ep = ridgePlanesHCP(ep);
        nAx = size(ep.ridgePlaneNormals,1);
    end
    for ii = 1:K
        results(ii,jj).eu     = eu(ii,:);
        results(ii,jj).label  = labels(jj);
        results(ii,jj).family = families{jj};
        results(ii,jj).drp    = DRPsim_hcp(eu(ii,1),eu(ii,2),eu(ii,3),ep);
        results(ii,jj).nAxes  = nAx;
    end
end

% ---- draw -------------------------------------------------------------
nCol = N + hasMeas;
fig = figure('Name','compareRidgeAxes', ...
    'Position',[100 100 min(1800,300*nCol) min(1000,260*K)]);
tl = tiledlayout(fig,K,nCol,'TileSpacing','compact','Padding','compact');

for ii = 1:K
    if hasMeas
        nexttile(tl,(ii-1)*nCol + 1);
        DRPdisp(DRP_norm(measured{ii}),exp_para);
        if ii == 1, title('measured','FontWeight','bold'); end
    end
    for jj = 1:N
        nexttile(tl,(ii-1)*nCol + hasMeas + jj);
        DRPdisp(results(ii,jj).drp,exp_para);
        if ii == 1
            title(sprintf('%s (%d)',labels(jj),results(ii,jj).nAxes), ...
                'Interpreter','none');
        end
    end
end
title(tl,sprintf('Ridge A/B  -  rows = orientations, mode = %s',mode), ...
    'FontWeight','bold');
end

% -------------------------------------------------------------------------
function lab = local_familyLabel(fam, mode)
if size(fam,1) == 1
    if mode == "direction"
        lab = sprintf('<%d %d %d %d>',fam(1),fam(2),fam(3),fam(4));
    else
        lab = sprintf('{%d %d %d %d}',fam(1),fam(2),fam(3),fam(4));
    end
else
    lab = sprintf('%d sets',size(fam,1));
end
end

% -------------------------------------------------------------------------
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
    error('compareRidgeAxes:measured', ...
        'measured must supply one DRP per orientation (K = %d).',K);
end
end
