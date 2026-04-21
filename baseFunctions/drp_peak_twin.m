% DRM Large-Scale Analysis

%% DRM measurement indexing engine
exp_para.th_max = 65;
exp_para.th_min = 15;
exp_para.th_num = 11;
exp_para.ph_num = 72;
exp_para.ph_min = 0;
exp_para.ph_max = 355;

exp_para.faceting = [1 0 0]; % once it gets rotated around, it should be fine

exp_para.fitting_para = [1 0 0 0];

pos1 = [1 1 2048 1080]; %[1 1 2048 1080]; %[580 1059 600 1079]; %one of the non_good groups takes


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

% turn the DRP into a summation for 3 and 4
drp_sum = sum(drp_4d,3);
drp_sum = squeeze(drp_sum);
[d1,d2,~] = size(drp_sum); 
map_peaks = zeros([d1,d2]);
peak_threshold = 3;
for i = 1:n1
    for j = 1:n2
        vec = squeeze(drp_sum(i,j,:));
        if max(vec) > 11*80
            [~, nr_peaks] = findpeaks(vec, MinPeakDistance=6, MinPeakHeight=0.6*11*80);
            % Process the number of peaks as needed
            if length(nr_peaks) > peak_threshold
                map_peaks(i,j) = 1;
            end
        end
     end     
end

figure;
imagesc(map_peaks);
colorbar;