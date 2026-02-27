%% EBSD stitching
rot = rotation.byAxisAngle(vector3d.X,0*degree);
ebsd_temp = rotate(ebsd,rot,'keepXY');
rot = rotation.byAxisAngle(vector3d.X,20*degree);
ebsd_temp = rotate(ebsd_20tilt,rot,'keepXY');
% ebsd_temp = ebsd_20tilt_r;
[grains, ebsd_temp.grainId] = calcGrains(ebsd_temp('indexed'),15*degree);
figure, plot(ebsd_temp('indexed'),ebsd_temp('indexed').orientations,'micronbar','off');
hold on
plot(grains.boundary,'LineWidth',2);
hold off
% original ebsd plot of the sample
%%
stepSize = max(ebsd_temp.unitCell) - min(ebsd_temp.unitCell);
pos_x = ebsd_temp.prop.x;
pos_y = ebsd_temp.prop.y;
pos_x = pos_x - min(pos_x);
pos_y = pos_y - min(pos_y);
idxEBSD = fix([pos_x, pos_y] ./ stepSize + [1,1]);
% under the current input settings, a 90-deg clockwise rotation is needed.
sizeEBSD = max(idxEBSD);
clear EUmap_temp
EUmap_temp(1,:,:) = rot90(reshape(ebsd_temp.rotations.phi1,sizeEBSD(1),sizeEBSD(2)),1)./degree;
EUmap_temp(2,:,:) = rot90(reshape(ebsd_temp.rotations.Phi,sizeEBSD(1),sizeEBSD(2)),1)./degree;
EUmap_temp(3,:,:) = rot90(reshape(ebsd_temp.rotations.phi2,sizeEBSD(1),sizeEBSD(2)),1)./degree;

isIndexedMap = rot90(reshape(ebsd_temp.isIndexed,sizeEBSD(1),sizeEBSD(2)),1);
grainIdMap = rot90(reshape(ebsd_temp.grainId,sizeEBSD(1),sizeEBSD(2)),1);

EUmap_temp(:,~isIndexedMap) = NaN;
EUmap_temp = permute(EUmap_temp,[2 3 1]);
color_ebsd_temp = plot_ipf_map(EUmap_temp,plotDir="z");
figure, imshow(color_ebsd_temp,'border','tight')
%% ----------------------------------------------------------------------------
% this part of transformation needs careful consideration
validRange = [floor(sizeEBSD(1)/4),floor(sizeEBSD(1)/4*3)];
shiftVal = (find(isIndexedMap(:,validRange(1)),1) - find(isIndexedMap(:,validRange(2)),1)) / (validRange(2)-validRange(1));
% shiftVal = 0;
tformMat = [1 0 0; shiftVal 1 0; 0 0 1];
tform = affinetform2d(tformMat);
EUmap_ebsd = imwarp(EUmap_temp,tform,'nearest');
isIndexedMap_ebsd = imwarp(isIndexedMap,tform,"nearest");
grainIdMap_ebsd = imwarp(grainIdMap,tform,"nearest");
EUmap_ebsd = permute(EUmap_ebsd,[3 1 2]);
EUmap_ebsd(:,~isIndexedMap_ebsd) = NaN;
EUmap_ebsd = permute(EUmap_ebsd,[2 3 1]);
figure, imshow(imwarp(color_ebsd_temp,tform,"nearest"),'Border','tight');
% ----------------------------------------------------------------------------
% plot scale bar according to EBSD dataset
% ----------------------------------------------------------------------------
%% if necessary, select the region of interest from plotted EBSD colormap
figure(101)
imshow(plot_ipf_map(EUmap_ebsd),'Border','tight')
roi_ebsd = drawrectangle;
pos_ebsd = roi_ebsd.Position;
close 101
EUmap_ebsd = imcrop(EUmap_ebsd,pos_ebsd);
figure, imshow(plot_ipf_map(EUmap_ebsd),'Border','tight')
isIndexedMap_ebsd = imcrop(isIndexedMap_ebsd,pos_ebsd);
grainIdMap_ebsd = imcrop(grainIdMap_ebsd,pos_ebsd);
%% register EBSD and DRM dataset and compare indexing error
% register two datasets
eumap = indexResult.euMap;
colorDRMoriginal = plot_ipf_map(eumap);
colorEBSDoriginal = plot_ipf_map(EUmap_ebsd);
if ~exist("movingPoints",'var')
    [movingPoints, refPoints] = cpselect(colorDRMoriginal,colorEBSDoriginal,'Wait',true);
end
tform_register = fitgeotrans(movingPoints,refPoints,'affine');
output_region = imref2d(size(colorEBSDoriginal));
EUmap_trans = imwarp(eumap,tform_register,'nearest','OutputView',output_region);
figure, imshowpair(plot_ipf_map(EUmap_trans),colorEBSDoriginal,'montage')

%% compare EUmap_trans and EUmap_ebsd
[n1,n2,~] = size(EUmap_ebsd);
eulerDRM = reshape(EUmap_trans,n1*n2,3);
eulerEBSD = reshape(EUmap_ebsd,n1*n2,3);
% cs = ebsd_temp.CSList{2};
cs = crystalSymmetry('cubic'); %change this
oriDRM = orientation.byEuler(eulerDRM.*degree,cs);
oriEBSD = orientation.byEuler(eulerEBSD.*degree,cs);
rot = rotation.byAxisAngle(vector3d.X,0*degree);
oriEBSD = rot*oriEBSD;
rot = rotation.byAxisAngle(vector3d.Y,30*degree);
oriEBSD = rot*oriEBSD;
rot = rotation.byAxisAngle(vector3d.Z,90*degree);
oriEBSD = rot*oriEBSD;

misOriAngle = angle(transpose(oriDRM),oriEBSD,cs)./degree;
misOriMap = abs(reshape(misOriAngle,n1,n2)-3);
misOriMap(~isIndexedMap_ebsd) = NaN;
% plot indexing error mapping
figure, imshow(misOriMap,Border="tight")
colormap(sky)
clim([0 20])
median(misOriMap(~isnan(misOriMap)),'all')
% colorbar
%% plot pixel-wise indexing error histogram
figure,
histogram(misOriMap,61,'BinLimits',[1 62],...
    'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
set(gca,'LineWidth',1.5,'FontSize',14)
xlabel('prediction error (deg)')
xlim([1 62])
ylim([0 5000])
ylabel('number of pixels')
title('histogram of orientation error')

%% calculate grain-wise error and plot corresponding histogram
grainSizeThres = 100;
num_grain = max(grainIdMap_ebsd,[],'all');
misOriGrain = zeros(num_grain,1);
for ii = 1:num_grain
    if sum(grainIdMap_ebsd == ii,"all") < grainSizeThres
        misOriGrain(ii) = NaN;
        continue
    else
        misOri_temp = misOriMap(grainIdMap_ebsd == ii);
        grainiiMap = grainIdMap_ebsd == ii;

        misOriGrain(ii) = median(misOri_temp(misOri_temp < prctile(misOri_temp,80)));
    end
end
figure,
histogram(misOriGrain,61,'BinLimits',[1 62],...
    'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
set(gca,'LineWidth',1.5,'FontSize',14)
% xlabel('prediction error (deg)')
xlim([1 62])
ylim([0 100])
% ylabel('number of grains')
% title('histogram of orientation error')
%%
figure, plotIPDF(oriEBSD,misOriAngle,vector3d.Z,'points',5e6,'MarkerSize',1)
colormap("jet")

%%
figure, plotIPDF([oriEBSD_top,oriEBSD_bot],[misOriAngle_top_02;misOriAngle_bot_02],vector3d.Z,'points',1e6,'MarkerSize',1)
colormap("jet")

%%
misOri = [reshape(misOriMap_top,[],1); reshape(misOriMap_bot,[],1)];
misOri = misOri(~isnan(misOri));
figure,
histogram(misOri,61,'BinLimits',[1 62],...
    'EdgeColor','k','EdgeAlpha',0.5,'FaceColor','#0072BD','FaceAlpha',1)
set(gca,'LineWidth',1.5,'FontSize',14)
xlabel('prediction error (deg)')
xlim([1 62])
ylim([0 1e5])
ylabel('number of pixels')
title('histogram of orientation error')

%% plot on indexing error IPF map
figure, plotIPDF(oriEBSD,misOriAngle,vector3d.Z,...
    'points',length(oriEBSD),'MarkerSize',1)
colormap("jet")


%% rotating back to sample coordinates
eulerDRM = reshape(EUmap_trans,n1*n2,3);
oriDRM = orientation.byEuler(eulerDRM.*degree,cs);
rot = rotation.byAxisAngle(vector3d.Y,-30*degree);
oriDRM = rot*oriDRM;
eulerDRM_sc = reshape([oriDRM.phi1, oriDRM.Phi, oriDRM.phi2]./degree, n1, n2, 3);
figure, imshow(plot_ipf_map(eulerDRM_sc),'Border','tight')
