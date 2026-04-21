function [drp_measurement, drp_predicted, x, y] = check_indexing_result(EUmap,drp_original,exp_para,options)
% this function is to compare DRM indexing result between measured DRP and
% predicted DRP from indexing result.
% create date: Sep 6, 2021
% edit date: Sep 20, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------

arguments
    EUmap double
    drp_original cell
    exp_para struct
    options.plot_xtal (1,1) logical = false
end

euler = reshape(EUmap,size(EUmap,1)*size(EUmap,2),size(EUmap,3));
% get the color mapping of DRM measurement
%cs = crystalSymmetry('cubic');%change this
%cs = crystalSymmetry('hex');
cs = {
    crystalSymmetry('6/mmm', [2.95, 2.95, 4.68], 'X||a*', 'Y||b', ...
    'mineral', 'Ti-Hex', 'color', 'light gray')
    };
oM = ipfHSVKey(cs);
ori_drm_all = orientation.byEuler(euler(:,1)*degree,euler(:,2)*degree,euler(:,3)*degree,cs);
color_drm_all = oM.orientation2color(ori_drm_all);
color_drm_reg_all = reshape(color_drm_all,size(EUmap,1),size(EUmap,2),3);
clear euler oM ori_drm_all color_drm_all

figure('Name','demo_fig');
imshow(color_drm_reg_all,'Border','tight');
%x = [52; 284; 18]; % for 400 400 411 411
%y = [110; 220; 284];
x = [90; 140; 170]; %Ti64
y = [40; 43; 14]; %Ti64
x = [];
y = [];
%fix it to get equal comparisons
%[x,y] = ginput;
% press 'enter' to stop
nn = length(y);
y = fix(y);
x = fix(x);
close(findobj('type','figure','name','demo_fig'));

%{
figure('Position',[200,200,200*(nn+2),200*2.5])
tiledlayout(2,nn+2,'TileSpacing','tight','Padding','compact')
nexttile(1,[2,2])
imshow(color_drm_reg_all,'Border','tight')
hold on
scatter(x,y,72,'x','k')
for ii = 1:nn
    text(x(ii)+5,y(ii)+5,int2str(ii),'FontSize',14)
end
hold off
%}
drpDic_q = evalin('base', 'drpDic');
euDic_q = evalin('base', 'euDic');

cs_q = cs{1};
oriLib = orientation.byEuler(euDic_q(:,1)*degree, euDic_q(:,2)*degree, euDic_q(:,3)*degree, cs_q);
figQuery = figure('Name', 'Euler Query');
axQuery = axes(figQuery);

fig = figure('Position',[200,200,400*(nn+2),400*2.5],'Name','check_indexing_result');
tl = tiledlayout(fig, 1, 3);
axMap = nexttile(tl,1);
imshow(color_drm_reg_all);
hold(axMap,'on');

axMeas = nexttile(tl,2);

axPred = nexttile(tl,3);

drp_measurement = cell(nn,1); %Ti64
drp_predicted = cell(nn,1); %Ti64
drp_measurement = {};
drp_predicted = {};
k = 0;
exitAll = false;

while isvalid(fig)
    figure(fig);
    axes(axMap);
    [x_pos, y_pos] = ginput(1);
    if isempty(x_pos) || isempty(y_pos)
        break
    end

    col = round(x_pos);
    row = round(y_pos);
    k = k+1;
    x(k,1) = col;
    y(k,1) = row;
    x_pos = fix(x_pos);
    y_pos = fix(y_pos);
    drp_measurement_tmp = drp_original{row,col};
    drp_measurement{k,1} = drp_measurement_tmp;

    cla(axMeas);
    axes(axMeas);
    DRPdisp(DRP_norm(drp_measurement_tmp),exp_para);

    eu_tmp = squeeze(EUmap(row,col,:)).';
    if exist('drpTable_1')
        if fitQuality(x_pos,y_pos) > length(drpTable_1)
            drpsim_tmp = DRPsim_double(eu_tmp(1),eu_tmp(2),eu_tmp(3),fitting_para(1:4));
        else
            drpsim_tmp = DRPsim_double(eu_tmp(1),eu_tmp(2),eu_tmp(3),[4,0.01,16,4]);
        end
    else
        drpsim_tmp = DRPsim(eu_tmp(1),eu_tmp(2),eu_tmp(3),exp_para);
    end

    cla(axPred);
    axes(axPred);
    DRPdisp(drpsim_tmp,exp_para);
    drp_predicted{k,1} = drpsim_tmp;
    eu_tmp

    drawnow;

    if isvalid(figQuery) && ~exitAll
        cla(axQuery);
        axes(axQuery);
        while true
            s = strtrim(input('Euler Query [phi1 Phi phi2] (Enter=continue, q=quit):','s'));
            if isempty(s)
                break
            end
            if strcmpi(s,'q')
                exitAll = true;
                break
            end

            vals = sscanf(strrep(s,',',''), '%f');
            if numel(vals) ~= 3
                fprintf('Please enter exactly 3 numbers');
                continue
            end

            q = vals(:).';
            oriQ = orientation.byEuler(q(1)*degree, q(2)*degree, q(3)*degree, cs_q);
            mis = angle(oriLib, oriQ, cs_q)./degree;
            [misMin, idxNearest] = min(mis);
            cla(axQuery); axes(axQuery);
            DRPdisp(drpDic_q{idxNearest}, exp_para);
            eNear = euDic_q(idxNearest,:);
            fprintf('Input: [%.2f %.2f %.2f]\n', q(1), q(2), q(3));
            fprintf('Nearest lib: [%.2f %.2f %.2f], id=%d, mis=%.3f deg\n', eNear(1), eNear(2), eNear(3), idxNearest, misMin);
        end
    end
end

if isvalid(fig)
    close(fig);
end

%{
for ii = 1:nn
    % DRP from measurement
    nexttile(ii+2)
    x_pos = y(ii);
    y_pos = x(ii);
    drp_measurement{ii} = drp_original{x_pos,y_pos};
    DRPdisp(DRP_norm(drp_measurement{ii}),exp_para)
    % DRP from prediction
    nexttile((ii+2+nn+2))
    eu_tmp = [EUmap(x_pos,y_pos,:)];
    if exist('drpTable_1')
        if fitQuality(x_pos,y_pos) > length(drpTable_1)
            drpsim_tmp = DRPsim_double(eu_tmp(1),eu_tmp(2),eu_tmp(3),fitting_para(1:4));
        else
            drpsim_tmp = DRPsim_double(eu_tmp(1),eu_tmp(2),eu_tmp(3),[4,0.01,16,4]);
        end
    else
        drpsim_tmp = DRPsim(eu_tmp(1),eu_tmp(2),eu_tmp(3),exp_para);
    end
    DRPdisp(drpsim_tmp,exp_para);
    drp_predicted{ii} = drpsim_tmp;
    eu_tmp
end
%}

if options.plot_xtal
    plot crystal shape
    figure('Position',[200,200,200*nn,200])
    tiledlayout(1,nn,'TileSpacing','tight','Padding','compact')
    for ii = 1:nn
        % DRP from measurement
        nexttile(ii)
        x_pos = y(ii);
        y_pos = x(ii);
        euler_angle = [EUmap(x_pos,y_pos,:)];
        %     cs = crystalSymmetry('cubic');
        %     ori_tmp = orientation.byEuler(eu_tmp(1)*degree, eu_tmp(2)*degree, eu_tmp(3)*degree,cs);
        %     cS = crystalShape.cube(cs);
        %     plot(ori_tmp * cS * 0.7);
        %     axis equal
        %     axis off
        eu1 = euler_angle(1);
        eu2 = euler_angle(2);
        eu3 = euler_angle(3);
        title(sprintf('%0.1f, %0.1f, %0.1f',[eu1, eu2, eu3]))

        drpsim_tmp = DRPsim(eu1,eu2,eu3,exp_para); % the matrix in form of th_num * ph_num
        C = [drpsim_tmp,drpsim_tmp(:,1)];
        %     theta = repmat([0:90/(exp_para.th_num-1):90]',1,exp_para.ph_num+1);
        %     phi = repmat([0:360/exp_para.ph_num:360],exp_para.th_num,1);
        %     x = cosd(theta).*cosd(phi);
        %     y = cosd(theta).*sind(phi);
        %     z = sind(theta);

        %     surf(x,y,z,C,'EdgeColor','none')
        %     axis equal
        %     set(gca,'visible','off')
        %     colormap('jet')
        %     alpha 0.45

        %change this faceting

        % Cube (for {100} faceting)
        %r = .3;
        %V1=[-1;  1; 1; -1; -1;  1; 1; -1;];
        %V2=[-1; -1; 1;  1; -1; -1; 1;  1;];
        %V3=[-1; -1;-1; -1;  1;  1; 1;  1;];
        %F= [1 2 3 4; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8;];
        %[THETA,PHI,R]=cart2sph(V1,V2,V3);
        %R=r.*ones(size(V1(:,1)));
        %[V1,V2,V3]=sph2cart(THETA,PHI,R);
        %V=[V1 V2 V3];
        % apply the rotation
        %roll = eu1*degree; pitch = eu2*degree; yaw = eu3*degree;
        %dcm = angle2dcm(roll, pitch, yaw, 'ZXZ');
        %V = V*dcm;
        %patch('Faces',F,'Vertices',V,'FaceColor',[77 77 77]/255,'FaceAlpha',0.6,'EdgeColor','k',...
        %    'LineWidth',2); axis equal; grid on; hold on; view(3);
        %axis off
        %view(0,20)

        % hexagonal for prismatic faceting
        %% Change this to the faceting that will be used
        r = .3;
        theta = (0:5)' * pi/3;

        V1 = [cos(theta); cos(theta)];
        V2 = [sin(theta); sin(theta)];
        V3 = [-ones(6,1); ones(6,1)];
        F = [1 2 7 8; 2 3 8 9; 3 4 9 10; 4 5 10 11; 5 6 11 12; 6 1 12 1];
        [THETA,PHI,R]=cart2sph(V1,V2,V3);
        R=r.*ones(size(V1(:,1)));
        [V1,V2,V3]=sph2cart(THETA,PHI,R);
        V=[V1 V2 V3];
        % apply the rotation
        roll = eu1*degree; pitch = eu2*degree; yaw = eu3*degree;
        dcm = angle2dcm(roll, pitch, yaw, 'ZXZ');
        V = V*dcm;
        patch('Faces',F,'Vertices',V,'FaceColor',[77 77 77]/255,'FaceAlpha',0.6,'EdgeColor','k',...
            'LineWidth',2); axis equal; grid on; hold on; view(3);
        axis off
        view(0,20)

    end
end