function [histogram_alignment] = drp_direction_analysis(drps)

% fold the drp to make phi1 and phi1+180 the same dimension
[n1, n2, th_num, ph_num] = size(drps);
size(drps)
folded_drps = drps(:,:,:,1:ph_num/2)+drps(:,:,:,ph_num/2+1:end);
sum_drps = sum(folded_drps,3);

[Cmax, IndC] = max(sum_drps,[],4);
max_phi_index = IndC;

phi_values = linspace(0,180,ph_num/2+1);
phi_values(end) = [];
max_phi_angle = phi_values(reshape(max_phi_index(:), [n1,n2]));

mean_peak_intensity = zeros(n1,n2);
for i = 1:n1
    for j = 1:n2
        peak_idx = max_phi_index(i,j);
        mean_peak_intensity(i,j) = mean(folded_drps(i,j,:,peak_idx),3);
    end
end

% plot the map
figDRM = figure('Name', 'DRM');
[n1,n2] = size(max_phi_angle)
max_phi_angle_90 = max_phi_angle - 90;
max_phi_angle_90(max_phi_angle_90 < 0) = ...
    max_phi_angle_90(max_phi_angle_90 < 0) + 180;
imagesc(max_phi_angle_90);
axis image;
colormap(jet);
hcb = colorbar;
hcb.Title.String = ["DRP Peak Direction (" char(176) ")"];
%title('Phi direction for maximum values');
title('DRP Peak Direction (°)')
axDRM = gca;
set(axDRM, 'YDir', 'reverse');

hold on;

figure;
imagesc(mean_peak_intensity);
axis image;
colormap(jet);
colorbar;
title('Mean Intensity Along Peak Direction')

figure;
histogram_alignment = histogram(max_phi_angle(:), phi_values);
title('Histogram distribution');

%choose landmarks to compare with the EBSD map
% turn max_phi_angle into ctf

trial_h = [74.19 82.36 81.21 78.14 85.40 89.92 82.60 87.93 87.33 77.13 86.09 18.08 20.63 78.64 8.89];
figure('Name', 'Trial Hisotgram');
histogram(trial_h);

%% Load EBSD data
ebsd = EBSD.load("C:\Users\mrbla\OneDrive\Bureaublad\Cambridge\P&W Deliverable 4\Ti7-Large-EBSD-Map-CTF.ctf");

% Restrict to Ti-Hex phase only
ebsd = ebsd('Ti-Hex');
[phi1, ~, ~] = Euler(ebsd.orientations);
figEBSD = figure('Name','EBSD');
plot(ebsd, phi1);
axEBSD = gca;
hold(axEBSD, 'on');

disp('Click one pixel in EBSD, then one in DRM. Close either window to stop.');

mislist = [];

%% Interactive loop
while isvalid(figEBSD) && isvalid(figDRM)
    % --- Select pixel in predicted EBSD ---
    f = figure(figDRM); drawnow;                   % ensure focus
    axes(axDRM);
    if f.CurrentCharacter > 0
        break;
    end
    [x1,y1, btn1] = ginput(1);
    if isempty(btn1) || btn1 ~=1, break; end
    col1 = round(x1);
    row1 = round(y1);
    phi1_DRM = max_phi_angle(row1,col1);
    scatter(axDRM, col1, row1, 60, 'g', 'filled');

    % --- Select pixel in true EBSD ---
    f = figure(figEBSD); drawnow;
    axes(axEBSD);
    if f.CurrentCharacter > 0
        break;
    end
    [x2,y2, btn2] = ginput(1);
    if isempty(btn2) || btn2 ~= 1, break; end
    [~,idx2] = min((ebsd.x - x2).^2 + (ebsd.y - y2).^2);
    ori_EBSD = ebsd.orientations(idx2);
    [e1, ~, ~] = Euler(ori_EBSD);
    phi1_EBSD = e1 * 180/pi; 
    scatter(axEBSD, ebsd.x(idx2), ebsd.y(idx2), 60, 'r', 'filled');

    % --- Misorientation ---
    mis = abs(phi1_DRM - phi1_EBSD);
    mis = min(mis, 180-mis);
    fprintf('Misorientation = %.2f°\n', mis);
    mislist = [mislist, mis];
end

figure('EBSD Histogram');
histogram(mislist);

end