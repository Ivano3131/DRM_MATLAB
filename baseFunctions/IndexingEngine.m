function result = IndexingEngine(drp_original,AE_DRM,exp_para,drpDic,euDic,rotDic)
% this is the indexing engine to index orientation of each pixels by its
% euler angle triplets. 
% The output are: (1) euler angle triplets mapping; (2)
% indexing quality (the method to be determined)
% Edit date: Sep 20, 2021
% By: Chenyang ZHU @ NTU
% -------------------------------------------------------------------------

% this relates the measured DRPs and the generated DRPs

[n1,n2] = size(drp_original);
[n3,n4] = size(drp_original{1});
th_num = exp_para.th_num;
ph_num = exp_para.ph_num;
if th_num ~= n3 || ph_num ~= n4
    error('size of drp stack does not match')
end
n_tot = n3 * n4;
num_dic = length(drpDic);
midpt = floor(ph_num/2);

dic_Enc = encode(AE_DRM,drpDic'); % encode the dictionary too

knn_num = 10;
% output
result.EUmap = zeros(n1,n2,3);
result.quality = zeros(n1,n2);
result.shift = zeros(n1,n2);
result.idx = zeros(n1,n2,knn_num);
% these four are pre-define output variable
for ii = 1:n1
    drp_row = drp_original(ii,:);
    shift_row = zeros(n2,1);
    drp_cell = cell(1,n2);
    parfor jj = 1:n2
        sum_col = sum(drp_row{jj} == max(max(drp_row{jj}))); % to be further considered
        [~, shift] = max(sum_col);
        shift = 0 - shift;
        tmpDRP = circshift(drp_row{jj},shift,2);
        shift_row(jj) = shift;
        drp_cell{jj} = double(tmpDRP)/255;
    end
    
    drp_Enc = encode(AE_DRM,drp_cell); % drp_Enc becomes a row of barcodes
    Idx = knnsearch(dic_Enc', drp_Enc','K',knn_num); %get the closest relation
    tmpEU = euDic(Idx(:,1),:);
    tmpEU(:,1) = (360 / ph_num) * (-shift_row + rotDic(Idx(:,1)));
    result.shift(ii,:) = shift_row;
    result.EUmap(ii,:,:) = tmpEU;
    result.quality(ii,:) = sqrt(sum((dic_Enc(:,Idx(:,1)) - drp_Enc).^2)); % difference between the Euler angles
    
    result.idx(ii,:,:) = Idx;
    workbar(ii/n1,sprintf('Computing texture mapping %d / %d',[ii n1]));
end

% drp_list = reshape(drp_original,n1*n2,1);
% shift_row = zeros(length(drp_list),1);
% parfor jj = 1:n1*n2
%     sum_col = sum(drp_list{jj});
%     [~, shift] = max(sum_col);
%     shift = midpt - shift;
%     tmpDRP = circshift(drp_list{jj},shift,2);
%     shift_row(jj) = shift;
%     drp_list{jj} = double(tmpDRP)/255;
% end
% % msgbox('calculation finished')
% drp_Enc = encode(AE_DRM,drp_list); % drp_Enc becomes a row of barcodes
% Idx = knnsearch(dic_Enc', drp_Enc');
% tmpEU = euDic(Idx,:);
% tmpEU(:,1) = (360 / ph_num) * (-shift_row + rotDic(Idx));
% % EUmap = zeros(n1,n2,3);
% % quality = zeros(n1,n2);
% EUmap = reshape(tmpEU,n1,n2,3);
% quality = sqrt(sum((dic_Enc(:,Idx) - drp_Enc).^2));
% quality = reshape(quality,n1,n2);
% msgbox('Indexing Completed!')

end