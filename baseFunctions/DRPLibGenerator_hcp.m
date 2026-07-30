function drpLib = DRPLibGenerator_hcp(ang_res, exp_para, options)
% DRPLibGenerator_hcp  DRP library for the direct (no-autoencoder) engine.
%
%   drpLib = DRPLibGenerator_hcp(ang_res, exp_para)
%
% Like makeDRPdic_hcp but keeps the DRPs unshifted (full DRP is matched
% directly by DirectDIEngine).  ang_res is in radians, e.g. 5*degree.
%
% Output
%   drpLib.drpDic    nn x 1 cell of th_num x ph_num uint8 DRPs
%   drpLib.eulerDic  nn x 3 Bunge Euler angles [phi1 Phi phi2] in degrees
% -------------------------------------------------------------------------
arguments
    ang_res (1,1) double
    exp_para struct
    options.verbose (1,1) logical = true
end

if ~isfield(exp_para,'facetNormals')
    exp_para = facetNormalsHCP(exp_para);
end
cs = exp_para.crystal;

ori = equispacedSO3Grid(cs,'resolution',ang_res);
nn  = length(ori);
phi1 = ori.phi1 ./ degree;
Phi  = ori.Phi  ./ degree;
phi2 = ori.phi2 ./ degree;

drpLib.drpDic   = cell(nn,1);
drpLib.eulerDic = zeros(nn,3);
for ii = 1:nn
    drpLib.drpDic{ii}     = DRPsim_hcp(phi1(ii),Phi(ii),phi2(ii),exp_para);
    drpLib.eulerDic(ii,:) = [phi1(ii) Phi(ii) phi2(ii)];
    if options.verbose
        workbar(ii/nn,sprintf("processing %d / %d DRPs",[ii nn]));
    end
end
end
