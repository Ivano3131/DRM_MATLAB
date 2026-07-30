function [all_rot, faceW, pairW] = rotate_facet_hcp(eu1,eu2,eu3,exp_para)
% rotate_facet_hcp  Specimen-frame HCP facet normals for one orientation.
%
%   [all_rot, faceW, pairW] = rotate_facet_hcp(eu1,eu2,eu3,exp_para)
%
% Rotates the precomputed CRYSTAL-frame facet normals (exp_para.facetNormals,
% produced by facetNormalsHCP) into the SPECIMEN frame for the Bunge Euler
% triplet (eu1,eu2,eu3) = (phi1,Phi,phi2) in degrees, and keeps the upward-
% pointing normals (the ones a facet on the top surface presents to the
% overhead camera).
%
% The rotation uses EulerRotate, i.e. specimen = Rz(phi1)*Rx(Phi)*Rz(phi2) *
% crystal, which is the same active Bunge convention as MTEX
% orientation.byEuler.  Because the crystal-frame normals come from the same
% MTEX crystalSymmetry object used for the orientation grid and the IPF map,
% the whole chain is on one consistent convention.
%
% This replaces the cubic rotate_facet.m: instead of permuting/sign-flipping a
% [u v w] vector (which is only valid in an orthonormal cubic frame), the
% crystallography is done once, correctly, in facetNormalsHCP.
%
% Outputs
%   all_rot  M x 3 unit facet normals in the specimen frame (z >= 0)
%   faceW    M x 1 peak weights   (copied from exp_para.faceW)
%   pairW    M x 1 band weights   (copied from exp_para.pairW)
% -------------------------------------------------------------------------
arguments
    eu1 double
    eu2 double
    eu3 double
    exp_para struct
end

if ~isfield(exp_para,'facetNormals') || isempty(exp_para.facetNormals)
    error('rotate_facet_hcp:notSetup', ...
        'exp_para.facetNormals is missing. Run exp_para = facetNormalsHCP(exp_para) first.');
end

n_xtal = exp_para.facetNormals;   % M x 3, crystal frame
faceW  = exp_para.faceW;
pairW  = exp_para.pairW;

% crystal -> specimen (Bunge, consistent with MTEX and the cubic pipeline)
all_rot = normr(EulerRotate(n_xtal, eu1, eu2, eu3));

% keep upward-pointing normals: a facet is only seen if its normal has a
% component towards the camera (+Z).
flipMask = all_rot(:,3) < 0;
all_rot(flipMask,:) = -all_rot(flipMask,:);
end
