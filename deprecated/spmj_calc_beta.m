function [beta,Sw]=spmj_calc_beta(Y,SPM,varargin);
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

params={};
M=[];
C0=[];
NP=1;

vararginoptions(varargin,{'params','C0','NP'});

% ---------------------------------------------------------------------
% Extract data and redo the GLM - get res'*res matrix
xX    = SPM.xX;                        %- Take the design
KWY   = spm_filter(xX.K,xX.W*Y);       %- Filter the data
clear Y;                               %- delete Y
beta  = xX.pKX*KWY;                    %- Parameter estimates
res   = KWY - xX.xKXs.X*beta;          %- Residuals - seems to be faster than spm_sp('r',xX.xKXs,KWY);
beta  = beta(xX.iC,:);                 %- get the task-relevant betas
T     = size(xX.X,1);                  
Sig   = res'*res./size(xX.X,1);          %- Sigma for all voxels
clear res KWY xX                       %- preserve memory

% Project out contrast of non-interest
if (~isempty(C0))
    beta  =  beta-C0*pinv(C0)*beta;
end;

[N,P]=size(beta);
if (P>T) 
    [U,S,V]=spm_svd(Sig); 
    iS=U*diag(sqrt(1./diag(S)))*V'; 
else
    iS1=Sig^(-0.5);                     % This is faster, but not secure 
end; 

betaS        = beta*iS;          % Find speed-up for this?
if (NP)
    result= feval(mva_func,betaS,params{:});
else
    result= feval(mva_func,betaS',params{:});
end;
