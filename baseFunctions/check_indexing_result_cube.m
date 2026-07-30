function [drp_measurement, drp_predicted, xy] = check_indexing_result_cube(EUmap, drp_original, exp_para, options)
% check_indexing_result_cube  Measured vs predicted DRPs, cube-lath model.
%
%   [drp_measurement, drp_predicted, xy] = ...
%       check_indexing_result_cube(EUmap, drp_original, exp_para)
%   ... = check_indexing_result_cube(EUmap, drp_original, exp_para, score=res.score)
%
% Shows the IPF-z map; click a pixel to see its measured DRP next to the DRP
% predicted by DRPsim_cube for the indexed orientation.  Press Enter (empty
% click) to finish.  Pass the match score to have it printed alongside - a
% low score with a plausible-looking pattern usually means the etch depth is
% wrong rather than the orientation.
%
% Outputs collect the clicked measured/predicted DRPs and their pixel
% coordinates xy = [col row].  Note drp_original is indexed {row,col} while
% the map reports (col,row).
%
% See also matchDRPcube, DRPsim_cube, checkCubeMechanism.
% -------------------------------------------------------------------------
arguments
    EUmap double
    drp_original cell
    exp_para struct
    options.score double = []
end

if ~isfield(exp_para,'cube') || isempty(exp_para.cube)
    error('check_indexing_result_cube:notSetup', ...
        'exp_para.cube is missing. Run exp_para = cubeGeometryTi64(exp_para) first.');
end

% IPF-z colour map of the indexing result
euler = reshape(EUmap, size(EUmap,1)*size(EUmap,2), size(EUmap,3));
cs = exp_para.crystal;
oM = ipfHSVKey(cs);
oM.inversePoleFigureDirection = vector3d.Z;
ori = orientation.byEuler(euler(:,1)*degree, euler(:,2)*degree, euler(:,3)*degree, cs);
rgb = reshape(oM.orientation2color(ori), size(EUmap,1), size(EUmap,2), 3);

fig = figure('Name','check_indexing_result_cube','Position',[150 150 1200 500]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
axMap  = nexttile(tl,1); imshow(rgb); title(axMap,'click a pixel (Enter to stop)');
axMeas = nexttile(tl,2); title(axMeas,'measured DRP');
axPred = nexttile(tl,3); title(axPred,'predicted DRP');
hold(axMap,'on');

drp_measurement = {};
drp_predicted   = {};
xy = zeros(0,2);
k = 0;

while isvalid(fig)
    [xc,yc] = ginput(1);
    if isempty(xc) || isempty(yc), break; end
    col = round(xc); row = round(yc);
    if row < 1 || col < 1 || row > size(EUmap,1) || col > size(EUmap,2)
        continue
    end
    k = k + 1;
    xy(k,:) = [col row];
    scatter(axMap,col,row,60,'k','x');

    meas = drp_original{row,col};
    eu   = squeeze(EUmap(row,col,:)).';
    [pred, E] = DRPsim_cube(eu(1),eu(2),eu(3),exp_para);

    cla(axMeas); axes(axMeas); DRPdisp(DRP_norm(meas),exp_para); title(axMeas,'measured DRP');
    cla(axPred); axes(axPred); DRPdisp(pred,exp_para);
    if isempty(options.score)
        title(axPred,sprintf('predicted  [%.1f %.1f %.1f]',eu(1),eu(2),eu(3)));
    else
        title(axPred,sprintf('predicted  [%.1f %.1f %.1f]   score %.3f', ...
            eu(1),eu(2),eu(3), options.score(row,col)));
    end

    % the numbers behind the picture: which edges the etch cut left standing
    fprintf('\npixel (col %d, row %d)  [phi1 Phi phi2] = [%.1f %.1f %.1f]\n', ...
        col, row, eu(1), eu(2), eu(3));
    fprintf('  %-12s %5s %7s %9s %8s %8s\n','edge','axis','tilt','stripeAz','expose','amp');
    for ii = 1:numel(E)
        fprintf('  %-12s %5d %7.1f %9.1f %8.3f %8.4f\n', ...
            E(ii).family, E(ii).axisId, E(ii).tilt, E(ii).stripeAz, ...
            E(ii).expose, E(ii).amp);
    end

    drp_measurement{k,1} = meas; %#ok<AGROW>
    drp_predicted{k,1}   = pred; %#ok<AGROW>
    drawnow;
end

if isvalid(fig), close(fig); end
end
