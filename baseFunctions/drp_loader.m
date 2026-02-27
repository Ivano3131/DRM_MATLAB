function [igrey,phitheta,pos,img_sample] = drp_loader(exp_para,pos,options)
% generate image stack variable as input for indexing engine.
% exp_para are global settings, compulsory
% ext_name is figure extension name, jpg as default
% scale_num is to scale up/down the original image, 1 as default
% pos_interest is ROI, optional
% Create date: Oct 27, 2022
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
    arguments
        exp_para struct
        pos (1,4) double
        options.format (1,:) string = 'jpg'
        options.scale (1,1) double = 1
        options.folder (1,:) string = ""
    end

    th_min = exp_para.th_min;
    th_max = exp_para.th_max;
    th_num = exp_para.th_num;
    ph_min = exp_para.ph_min;
    ph_max = exp_para.ph_max;
    ph_num = exp_para.ph_num;

    if options.folder == ""
        image_folder = uigetdir('*.*','Select folder containing images');
    else
        image_folder = options.folder;
    end

    ext = strcat('*.', options.format);
    D_im = dir(fullfile(image_folder, ext));

    num_img = length(D_im);
    if num_img ~= th_num * ph_num
        error('experimental parameters are not correct!')
    end

    scale = options.scale;
    tmp_file = fullfile(D_im(num_img).folder,D_im(num_img).name);
    tmp_img = imread(tmp_file);
    if length(size(tmp_img)) == 3
        tmp_img = imresize(rgb2gray(imread(tmp_file)),scale);
    else
        tmp_img = imresize(imread(tmp_file),scale);
    end
    % crop the image based on roi position

    if pos == zeros(1,4)
        figure(101)
        imshow(tmp_img * 2,'Border','tight')
        roi = drawrectangle;
        pos = roi.Position;
        close 101
    end

    % generate the phitheta angle profile
    % phi from ph_min -> ph_max for each theta row
    % theta from th_min -> th_max
    phitheta = zeros(num_img,2);
    indexing = (1:num_img)';
    phi_step = (ph_max - ph_min) / (ph_num - 1);
    th_step = (th_max - th_min) / (th_num - 1);

    for ii = 1:num_img
        phitheta(ii,1) = rem(indexing(ii)-1,ph_num) * phi_step;
        phitheta(ii,2) = floor((ii-1)/ph_num)*th_step + th_min;
    end

    tmp_img_crop = imcrop(tmp_img,pos);
    img_sample = tmp_img_crop;
    [n1,n2] = size(tmp_img_crop);
    igrey = zeros(n1,n2,num_img,'uint8');
    if length(size(imread(tmp_file))) == 3
        for ii = 1:num_img
            tmp_file = fullfile(D_im(ii).folder,D_im(ii).name);
            igrey_tmp = imresize(rgb2gray(imread(tmp_file)),scale);
            name = split(D_im(ii).name,["_","."]);
            phi_tmp = str2num(name{1});
            theta_tmp = str2num(name{2})/10;
            %phitheta
            %phi_tmp
            %theta_tmp
            idx = find(phitheta(:,1)==phi_tmp & phitheta(:,2)==theta_tmp);
            %length(idx)
            igrey(:,:,idx) = imcrop(igrey_tmp,pos);
            workbar(ii/num_img,sprintf('Generating igrey, %d / %d',[ii num_img]));
        end
    else
        for ii = 1:num_img
            tmp_file = fullfile(D_im(ii).folder,D_im(ii).name);
            igrey_tmp = imresize(imread(tmp_file),scale);
            name = split(D_im(ii).name,["_","."]);
            phi_tmp = str2num(name{1});
            theta_tmp = str2num(name{2}); % /10 was here
            %phitheta
            %phi_tmp
            %theta_tmp
            idx = find(phitheta(:,1)==phi_tmp & phitheta(:,2)==theta_tmp);
            %length(idx)
            igrey(:,:,idx) = imcrop(igrey_tmp,pos);
            workbar(ii/num_img,sprintf('Generating igrey, %d / %d',[ii num_img]));
        end
    end

end