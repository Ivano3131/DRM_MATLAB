function [drp_measurement, drp_predicted, xy] = check_indexing_result_Ti64(EUmap,drp_original,exp_para)
% check_indexing_result_Ti64  Measured vs predicted DRP at clicked pixels.
%
%   [drp_measurement, drp_predicted, xy] = ...
%       check_indexing_result_Ti64(EUmap, drp_original, exp_para)
%
% Same interaction as check_indexing_result_hcp, but the prediction comes from
% DRPsim_Ti64 (the alpha-lath etching mechanism).  Shows the IPF-z map of the
% indexing result; click a pixel to see its measured DRP next to the predicted
% one, and the predicted stripe azimuths are printed to the command window so
% the mechanism can be checked pixel by pixel.  Press Enter (empty click) to
% finish.
%
% Outputs collect the clicked measured/predicted DRPs and the pixel
% coordinates xy = [col row].
% -------------------------------------------------------------------------
arguments
    EUmap double
    drp_original cell
    exp_para struct
end

if ~isfield(exp_para,'lath') || isempty(exp_para.lath)
    error('check_indexing_result_Ti64:notSetup', ...
        'exp_para.lath is missing. Run exp_para = lathGeometryTi64(exp_para) first.');
end

% IPF-z colour map of the indexing result
euler = reshape(EUmap,size(EUmap,1)*size(EUmap,2),size(EUmap,3));
cs = exp_para.crystal;
oM = ipfHSVKey(cs);
oM.inversePoleFigureDirection = vector3d.Z;
ori = orientation.byEuler(euler(:,1)*degree,euler(:,2)*degree,euler(:,3)*degree,cs);
rgb = reshape(oM.orientation2color(ori),size(EUmap,1),size(EUmap,2),3);

fig = figure('Name','check_indexing_result_Ti64','Position',[150 150 1200 500]);
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
    [pred, R] = DRPsim_Ti64(eu(1),eu(2),eu(3),exp_para);

    cla(axMeas); axes(axMeas); DRPdisp(DRP_norm(meas),exp_para); title(axMeas,'measured DRP');
    cla(axPred); axes(axPred); DRPdisp(pred,exp_para);
    title(axPred,sprintf('predicted  [%.1f %.1f %.1f]',eu(1),eu(2),eu(3)));

    % report the mechanism behind the prediction
    fprintf('pixel (row %d, col %d)  [phi1 Phi phi2] = [%.1f %.1f %.1f]\n', ...
        row, col, eu(1), eu(2), eu(3));
    for ii = 1:numel(R)
        f = R(ii);
        if f.kind == "ridge"
            fprintf('   %-14s stripe az %6.1f deg   amp %.4f\n', f.name, f.stripeAz, f.amp);
        else
            fprintf('   %-14s peak theta %6.1f deg  amp %.4f\n', f.name, f.peakTh, f.amp);
        end
    end

    drp_measurement{k,1} = meas; %#ok<AGROW>
    drp_predicted{k,1}   = pred; %#ok<AGROW>
    drawnow;
end

if isvalid(fig), close(fig); end
end
