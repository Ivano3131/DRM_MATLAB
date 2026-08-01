function [drp_uni] = splitDRP(exp_para,long_drp,phitheta)
% split a 1-d drp array into normal 2-d drp matrix
% Edit date: Aug 27, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------
ph_num = exp_para.ph_num;
th_num = exp_para.th_num;
% 'like' rather than a hard-coded uint8: a .drp carries uint16, and filling a
% uint8 buffer from it would saturate every sample above 255 to white.  For a
% uint8 stack this is the same array as before.
drp = zeros(ph_num,th_num,'like',long_drp);
th_step = (exp_para.th_max - exp_para.th_min) / (th_num - 1);
ph_step = 360 / ph_num;

for index = 1:length(phitheta)
    x = phitheta(index,1) / ph_step + 1;
    y = (phitheta(index,2) - exp_para.th_min) / th_step + 1;
    drp(x,y) = long_drp(index);
end

drp_uni = transpose(drp);

end