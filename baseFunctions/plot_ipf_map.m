function color_drm_reg = plot_ipf_map(EUmap,options)

arguments
    EUmap
    options.plotDir (1,1) string = "z"
end
euler = reshape(EUmap,size(EUmap,1)*size(EUmap,2),size(EUmap,3));
% get the color mapping of DRM measurement
cs = crystalSymmetry('cubic'); % change this
oM = ipfHSVKey(cs);
ori_drm = orientation.byEuler(euler(:,1)*degree,euler(:,2)*degree,euler(:,3)*degree,cs);
if options.plotDir == "z"
    oM.inversePoleFigureDirection = vector3d.Z;
elseif options.plotDir == "y"
    oM.inversePoleFigureDirection = vector3d.Y;
end
color_drm = oM.orientation2color(ori_drm);
color_drm_reg = reshape(color_drm,size(EUmap,1),size(EUmap,2),3);
% figure, imshow(color_drm_reg,'Border','tight')

end