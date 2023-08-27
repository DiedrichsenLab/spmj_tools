function [beta,Sw]=spmj_calc_beta_Sw(Y,xX,partition);
% function D=lmva_roi_raw(mva_func,SPM,ROI,varargin);
% Runs a mva-method on the beta's of a glm after using the first-level
% residuals to spatially prewhiten the data.
% INPUT:
%   Y:              
%   SPM:           SPM Structure
% VARARGIN:
%   params:        extra parameters passed to the mva-function
% OUTPUT:
%   result:        result vector (usually column) of the MVA function
% JD: voxel assignment for smaller number of subspaces guaranteed


% ---------------------------------------------------------------------
% Extract data and redo the GLM - get res'*res matrix
KWY   = spm_filter(xX.K,xX.W*Y);       %- Filter the data
clear Y;                               %- delete Y
beta  = xX.pKX*KWY;                    %- Parameter estimates
res   = KWY - xX.xKXs.X*beta;          %- Residuals - seems to be faster than spm_sp('r',xX.xKXs,KWY);

Sig   = res'*res./size(xX.X,1);          %- Sigma for all voxels
clear KWY xX                       %- preserve memory

part=unique(partition)'; 
for i=1:length(part); 
    indx=find(partition==part(i));
    Sw(:,:,i)=res(indx,:)'*res(indx,:); 
    Sw(:,:,i)=Sw(:,:,i)/length(indx);
end; 
