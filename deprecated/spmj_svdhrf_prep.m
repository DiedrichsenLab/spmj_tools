function [SPM,SPMconv]=spmj_svdhrf_prep(SPM,name,length,order,duration);
% function [SPMfir,xXconv]=spmj_svdhrf_prep(SPM);
% Takes a SPM structure and prepares it for
% HRF estimation (brain, region, or voxel-wise)
% INPUT:
%       SPM: normal SPM structure
%       name: What type of basis function set 
%       length: Length is s 
%       order: degre 
%       duration: duration of the events (optional) 
% OUTPUT:
%       SPMfir: copy of the SPM-structure with a finite response
%               function basis
%       xXconv: copy of the xX structure, prepared for convolution with the
%               estimate hrf
% Joern Diedrichsen


% Make an finite impulse response version of the old design
% matrix
SPM.xBF.name=name;
SPM.xBF.length=length;
SPM.xBF.order=order;
SPM.xBF=spm_get_bf(SPM.xBF);
if (nargin<5) 
    SPM=spmj_fMRI_design_changeBF(SPM,0);
else 
    SPM=spmj_fMRI_design_changeBF(SPM,duration);
end; 
xX    = SPM.xX;                        %- Take the design

% Extract the convolution matrix
i         = find(xX.iCb==1 | xX.iN>0);
xX.X   = xX.X(:,i);
xX.Xv  = xX.X; % This is the copy that will be changed on a voxel-by-voxel basis
xX.iCc = xX.iCc(i);
xX.iN  = xX.iN(i);
xX.iCs = xX.iCs(i);
xX.iC  = find(xX.iCc>0 | xX.iN==1);
xX.iB  = find(xX.iN==2);                 % constant terms
SPMconv    = SPM; 
SPMconv.xX = xX; 