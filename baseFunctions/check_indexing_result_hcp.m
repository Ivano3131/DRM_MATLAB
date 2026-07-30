function [drp_measurement, drp_predicted, xy] = check_indexing_result_hcp(EUmap,drp_original,exp_para)
% check_indexing_result_hcp  Compare measured vs predicted DRPs (HCP).
%
%   [drp_measurement, drp_predicted, xy] = ...
%       check_indexing_result_hcp(EUmap, drp_original, exp_para)
%
% Shows the IPF-z map; click a pixel to see its measured DRP next to the DRP
% predicted by DRPsim_hcp for the indexed orientation.  Press Enter (empty
% click) to finish.  This is the primary tool for checking, on a grain of
% known orientation, whether the chosen facet family reproduces the real
% peak/band geometry.
%
% Outputs collect the clicked measured/predicted DRPs and their pixel
% coordinates xy = [col row].
% -------------------------------------------------------------------------
arguments
    EUmap double
    drp_original cell
    exp_para struct
end

if ~isfield(exp_para,'facetNormals')
    exp_para = facetNormalsHCP(exp_para);
end

% IPF-z colour map of the indexing result
euler = reshape(EUmap,size(EUmap,1)*size(EUmap,2),size(EUmap,3));
cs = exp_para.crystal;
oM = ipfHSVKey(cs);
oM.inversePoleFigureDirection = vector3d.Z;
ori = orientation.byEuler(euler(:,1)*degree,euler(:,2)*degree,euler(:,3)*degree,cs);
rgb = reshape(oM.orientation2color(ori),size(EUmap,1),size(EUmap,2),3);

fig = figure('Name','check_indexing_result_hcp','Position',[150 150 1200 500]);
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
    pred = DRPsim_hcp(eu(1),eu(2),eu(3),exp_para);

    cla(axMeas); axes(axMeas); DRPdisp(DRP_norm(meas),exp_para); title(axMeas,'measured DRP');
    cla(axPred); axes(axPred); DRPdisp(pred,exp_para);
    title(axPred,sprintf('predicted  [%.1f %.1f %.1f]',eu(1),eu(2),eu(3)));

    drp_measurement{k,1} = meas; %#ok<AGROW>
    drp_predicted{k,1}   = pred; %#ok<AGROW>
    drawnow;
end

if isvalid(fig), close(fig); end
end
