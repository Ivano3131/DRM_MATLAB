function results = compareFacetFamilies(eu, exp_para, options)
% compareFacetFamilies  A/B simulated DRPs across candidate facet families.
%
%   results = compareFacetFamilies(eu, exp_para)
%   results = compareFacetFamilies(eu, exp_para, families=..., labels=..., ...
%                                  measured=..., peakModel=...)
%
% For one or more KNOWN orientations, simulates the DRP with DRPsim_hcp for
% each candidate {hkil} facet family and tiles them side by side, so the family
% that best reproduces a real (measured) DRP can be picked by eye.  This is the
% recommended way to decide exp_para.facetHKIL for a new material/etchant.
%
% Inputs
%   eu        K x 3 Bunge Euler angles [phi1 Phi phi2] in degrees (one per row).
%             Get these from EBSD on the same grains you have measured DRPs for.
%   exp_para  struct with the usual DRP grid + fitting_para fields and (ideally)
%             exp_para.crystal.  facetHKIL is overridden per candidate.
%
% Options
%   families  1 x N cell; each cell an M x 4 {h k i l} matrix (a candidate family
%             or a combined set).  Default: the four common alpha-Ti candidates.
%   labels    1 x N string; column titles.  Default: derived from families.
%   measured  K x 1 cell (or th_num x ph_num x K array) of measured DRPs to show
%             in a leftmost reference column.  Optional.
%   peakModel "specular" (default) or "cosine"; applied to all candidates.
%
% Output
%   results   struct array (K x N) with fields .eu .label .hkil .drp .nFacets
% -------------------------------------------------------------------------
arguments
    eu double
    exp_para struct
    options.families cell = { [1 0 -1 0], [1 1 -2 0], [1 0 -1 1], [0 0 0 1] }
    options.labels string = string([])
    options.measured = []
    options.peakModel string = "specular"
end

if size(eu,2) ~= 3
    error('compareFacetFamilies:eu','eu must be K x 3 [phi1 Phi phi2] in degrees.');
end
K = size(eu,1);
families = options.families;
N = numel(families);

% default column labels from the {hkil} of each candidate
if isempty(options.labels)
    labels = strings(1,N);
    for jj = 1:N
        labels(jj) = local_familyLabel(families{jj});
    end
else
    labels = options.labels;
end

% normalise the optional measured DRPs into a K x 1 cell
measured = local_measuredToCell(options.measured, K);
hasMeas = ~isempty(measured);

exp_para.peakModel = options.peakModel;

11

% pre-simulate everything (so the figure draws in one pass)
results = struct('eu',[],'label',[],'hkil',[],'drp',[],'nFacets',[]);
for jj = 1:N
    ep = exp_para;
    ep.facetHKIL = families{jj};
    ep.facetFaceWeights = [];   % uniform weights for a fair comparison
    ep.facetBandWeights = [];
    ep = facetNormalsHCP(ep);

    12
    
    for ii = 1:K
        results(ii,jj).eu      = eu(ii,:);
        results(ii,jj).label   = labels(jj);
        results(ii,jj).hkil    = families{jj};
        results(ii,jj).drp     = DRPsim_hcp(eu(ii,1),eu(ii,2),eu(ii,3),ep);
        results(ii,jj).nFacets = size(ep.facetNormals,1);
    end
end

% ---- draw -------------------------------------------------------------
nCol = N + hasMeas;
fig = figure('Name','compareFacetFamilies', ...
    'Position',[100 100 min(1800,300*nCol) min(1000,260*K)]);
tl = tiledlayout(fig,K,nCol,'TileSpacing','compact','Padding','compact');

for ii = 1:K
    if hasMeas
        nexttile(tl,(ii-1)*nCol + 1);
        DRPdisp(DRP_norm(measured{ii}),exp_para);
        if ii == 1, title('measured','FontWeight','bold'); end
        ylabel(sprintf('[%.0f %.0f %.0f]',eu(ii,:)));  % note: axis is off in DRPdisp
    end
    for jj = 1:N
        nexttile(tl,(ii-1)*nCol + hasMeas + jj);
        DRPdisp(results(ii,jj).drp,exp_para);
        if ii == 1
            title(sprintf('%s (%d)',labels(jj),results(ii,jj).nFacets), ...
                'Interpreter','none');
        end
    end
end
title(tl,sprintf('Facet-family A/B  -  rows = orientations, peakModel = %s', ...
    options.peakModel),'FontWeight','bold');
end

% -------------------------------------------------------------------------
function lab = local_familyLabel(hkil)
% short label like "{1 0 -1 0}" for a single family, or "2 families" for a set
if size(hkil,1) == 1
    lab = sprintf('{%d %d %d %d}',hkil(1),hkil(2),hkil(3),hkil(4));
else
    lab = sprintf('%d families',size(hkil,1));
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
else                        % single 2-D DRP for K==1
    measured = {m};
end
if numel(measured) ~= K
    error('compareFacetFamilies:measured', ...
        'measured must supply one DRP per orientation (K = %d).',K);
end
end
