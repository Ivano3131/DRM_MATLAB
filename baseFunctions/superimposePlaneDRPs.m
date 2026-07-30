function out = superimposePlaneDRPs(eu, exp_para, planes, options)
% superimposePlaneDRPs  Decompose/superimpose the DRP contributions of several
%                       crystallographic planes for a given orientation.
%
%   out = superimposePlaneDRPs(eu, exp_para, planes)
%   out = superimposePlaneDRPs(eu, exp_para, planes, measured=..., combine=...)
%
% Tests the multi-plane hypothesis: that a measured DRP is the superposition of
% stripes/dots from more than one crystallographic plane (e.g. basal + one
% selected prism/pyramidal variant), with orientation-dependent visibility.
% For each plane it simulates that plane's own DRP with DRPsim_hcp, then shows
% every contribution separately alongside the combined pattern and (optionally)
% the measured DRP.
%
% Inputs
%   eu       K x 3 Bunge Euler angles [phi1 Phi phi2] in degrees (rows = grains).
%   exp_para struct (uses .crystal, DRP grid, fitting fields).
%   planes   1 x N cell.  Each entry is either a 1x4 {h k i l} (defaults to a
%            ridge/family contribution, weight 1) OR a cell/struct specifying:
%              {hkil, model, variant, weight}
%                hkil    1 x 4 {h k i l} (Miller-Bravais plane)
%                model   "ridge" (surface-trace stripe, default) | "facet" (dot)
%                variant "family" (all symmetry equivalents, default) |
%                        "single" (that exact plane only -> variant selection)
%                weight  relative brightness in the combined pattern (default 1)
%
% Options
%   measured  K x 1 cell (or th x ph x K) of measured DRPs (optional).
%   combine   "max" (default, matches how DRPsim builds a DRP) | "sum"
%             (incoherent addition of independent features).
%
% Output
%   out  struct with fields:
%          .contrib   K x N cell of per-plane DRPs (double, [0,1])
%          .combined  K x 1 cell of the combined DRP (double, [0,1])
%          .labels    1 x N string labels
% -------------------------------------------------------------------------
arguments
    eu double
    exp_para struct
    planes cell
    options.measured = []
    options.combine string = "max"
end

K = size(eu,1);
N = numel(planes);
cs = exp_para.crystal;

if isempty(options.measured)
    meas = {};
else
    meas = local_measuredToCell(options.measured, K);
end
hasMeas = ~isempty(meas);

% ---- build one exp_para template per plane ----------------------------
epList  = cell(1,N);
weights = ones(1,N);
labels  = strings(1,N);
for jj = 1:N
    [hkil, model, variant, w] = local_parsePlane(planes{jj});
    weights(jj) = w;
    labels(jj)  = sprintf('{%d %d %d %d} %s/%s',hkil(1),hkil(2),hkil(3),hkil(4),model,variant);

    ep = exp_para;
    ep.reflectorModel = model;
    if model == "facet"
        if variant == "single"
            ep.facetNormals = local_singleNormal(hkil,cs);
            ep.faceW = 1; ep.pairW = 1;
        else
            ep.facetHKIL = hkil; ep.facetFaceWeights = []; ep.facetBandWeights = [];
            ep = facetNormalsHCP(ep);
        end
    else   % ridge (planeTrace)
        ep.ridgeMode = "planeTrace";
        if variant == "single"
            ep.ridgePlaneNormals = local_singleNormal(hkil,cs);
            ep.ridgeW = 1;
        else
            ep.ridgeHKIL = hkil; ep.ridgeWeights = [];
            ep = ridgePlanesHCP(ep);
        end
    end
    epList{jj} = ep;
end

% ---- simulate + combine ----------------------------------------------
out.contrib  = cell(K,N);
out.combined = cell(K,1);
out.labels   = labels;
for g = 1:K
    combined = [];
    for jj = 1:N
        drp = double(DRPsim_hcp(eu(g,1),eu(g,2),eu(g,3),epList{jj})) / 255;
        out.contrib{g,jj} = drp;
        c = weights(jj) * drp;
        if isempty(combined)
            combined = c;
        elseif options.combine == "sum"
            combined = combined + c;
        else
            combined = max(combined, c);
        end
    end
    if ~isempty(combined) && max(combined(:)) > 0
        combined = combined / max(combined(:));
    end
    out.combined{g} = combined;
end

% ---- draw -------------------------------------------------------------
nCol = hasMeas + N + 1;                     % [measured] contributions [combined]
fig = figure('Name','superimposePlaneDRPs', ...
    'Position',[80 80 min(1900,240*nCol) min(1050,240*K)]); %#ok<NASGU>
tl = tiledlayout(K,nCol,'TileSpacing','compact','Padding','compact');
for g = 1:K
    col = 0;
    if hasMeas
        col = col + 1; nexttile(tl,(g-1)*nCol+col);
        DRPdisp(DRP_norm(meas{g}),exp_para);
        if g==1, title('measured','FontWeight','bold'); end
    end
    for jj = 1:N
        col = col + 1; nexttile(tl,(g-1)*nCol+col);
        DRPdisp(out.contrib{g,jj},exp_para);
        if g==1, title(labels(jj),'Interpreter','none'); end
    end
    col = col + 1; nexttile(tl,(g-1)*nCol+col);
    DRPdisp(out.combined{g},exp_para);
    if g==1, title(sprintf('combined (%s)',options.combine),'FontWeight','bold'); end
end
title(tl,'Multi-plane superposition  -  rows = orientations','FontWeight','bold');
end

% -------------------------------------------------------------------------
function [hkil, model, variant, w] = local_parsePlane(p)
model = "ridge"; variant = "family"; w = 1;
if isnumeric(p)
    hkil = p;
elseif iscell(p)
    hkil = p{1};
    if numel(p) >= 2 && ~isempty(p{2}), model   = string(p{2}); end
    if numel(p) >= 3 && ~isempty(p{3}), variant = string(p{3}); end
    if numel(p) >= 4 && ~isempty(p{4}), w       = p{4};         end
elseif isstruct(p)
    hkil = p.hkil;
    if isfield(p,'model')   && ~isempty(p.model),   model   = string(p.model);   end
    if isfield(p,'variant') && ~isempty(p.variant), variant = string(p.variant); end
    if isfield(p,'weight')  && ~isempty(p.weight),  w       = p.weight;          end
else
    error('superimposePlaneDRPs:badPlane','each plane must be numeric, cell or struct.');
end
if numel(hkil) ~= 4
    error('superimposePlaneDRPs:badHKIL','each plane hkil must be 1 x 4.');
end
hkil = hkil(:).';
end

% -------------------------------------------------------------------------
function v = local_singleNormal(hkil, cs)
m = Miller(hkil(1),hkil(2),hkil(3),hkil(4),cs,'hkl');
v = normr(m.xyz);      % single crystal-frame plane normal, no symmetrise
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
    error('superimposePlaneDRPs:measured', ...
        'measured must supply one DRP per orientation (K = %d).',K);
end
end
