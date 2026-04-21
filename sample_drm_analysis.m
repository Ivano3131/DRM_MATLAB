% DRM Large-Scale Analysis

%% DRM measurement indexing engine
exp_para.th_max = 60;
exp_para.th_min = 20;
exp_para.th_num = 5;
exp_para.ph_num = 72;
exp_para.ph_min = 0;
exp_para.ph_max = 355;

exp_para.faceting = [1 0 0]; % once it gets rotated around, it should be fine

exp_para.fitting_para = [1 0 0 0];

pos1 = [500 200 1900 1000]; %[1 1 2048 1080]; %[580 1059 600 1079]; %one of the non_good groups takes


%% load sample and background dataset
scaleCoeff = 0.5;

%load the pictures dataset
[igray_sample, phitheta, pos, img_sample] = drp_loader( ...
    exp_para,pos1,format='jpg',scale=scaleCoeff);

%convert to DRP
drp_original = igray2drp(igray_sample,phitheta,exp_para); %[n1xn2xth_numxph_num]
[n1,n2] = size(drp_original);
drp_4d = zeros(n1,n2,exp_para.th_num,exp_para.ph_num, 'like', drp_original{1,1});
for i = 1:n1
    for j = 1:n2
        drp_4d(i,j,:,:) = drp_original{i,j};
    end
end
%drp_original = cell2mat(drp_original);

%get the average intensity value for each tile
intensity_map = squeeze(mean(mean(drp_4d,4),3));

%plot it as a histogram
intensity_distribution = reshape(intensity_map, [],1);
% Create a histogram of the intensity distribution
figure;
histogram(intensity_distribution, 'Normalization', 'probability');
xlabel('Intensity Value');
ylabel('Probability');
title('Intensity Distribution Histogram');

%maximum
intensity_max_map = squeeze(max(max(drp_4d, [], 4), [],3));
max_distribution = intensity_max_map(:);

figure;
histogram(max_distribution, 'Normalization', 'probability')
xlabel('Maximum intensity value');
ylabel('Probability');
title('Maximum Intensity Distribution');

%plot it as a map on top of one of the figures
figure;
imagesc(img_sample);
axis image;
colorbar;
title('Mean Intensity Distribution')

figure;
imagesc(intensity_max_map);
axis image;
colormap(jet);
colorbar;
title('Max Intensity Map')

% plot the maximum intensity distribution divided into five bins
numBins = 6;
% Create bins for the maximum intensity distribution
edges = [0 50 100 150 200 250 257];

binned_intensity_map = discretize(intensity_max_map, edges);
figure;
imagesc(binned_intensity_map);
axis image;
colorbar;
title('Binned Intensity Map')

% create the phi_1 - DRP direction map
histogram_phi_euler = drp_direction_analysis(drp_4d);