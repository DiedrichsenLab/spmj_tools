function [beta,ResSS,y_hat,y_res,r2]=spmj_svdhrf_est_roi(Y,SPM,SPMconv);
% function [SPMfir,xXconv]=spmj_svdhrf_prep(SPM);
% Takes a SPM structure and prepares it for
% HRF estimation (brain, region, or voxel-wise)
% INPUT:
%       Y:      N x P matrix of raw data 
%       SPM:     Design structure for FIR estimation 
%       xXconv: Design structure for convolution (both from spmj_svdhrf_prep)
% OUTPUT:
%       beta:   Estimated beta (QxP) 
%       ResSS:  residual sum of square (1xP)  
%       HRF:    the HRF that was used  (numFIRx1)
%       res:    full NxP residuals
% Joern Diedrichsen
xX    = SPM.xX; 
KWY   = spm_filter(xX.K,xX.W*Y);       %- Filter the data
beta  = xX.pKX*KWY;                    %- Parameter estimates
numbasis = max(xX.iCb);
P     = size(beta,2);
B=reshape(beta(xX.iCb>0,:),numbasis,sum(xX.iCb==1)*P);
[v,a]=svds(B,1);
SPMconv.xBF.bf=SPM.xBF.bf*v;
SPMconv.xBF.bf=SPMconv.xBF.bf./sum(SPMconv.xBF.bf); 
SPMconv.xBF.name='fitted';
SPMconv=spmj_fMRI_design_changeBF(SPMconv); 

% Now estimate the betas based on the estimated HRF
xXconv         = SPMconv.xX;
beta           = xXconv.pKX*KWY;
y_res          = spm_sp('r',xXconv.xKXs,KWY);       %-Residuals
ResSS          = sum(y_res.^2);                   %-Residual SS
r2             = 1-sum(KWY.^2)./ResSS; 

reg_interest=[xXconv.iH xXconv.iC]; 
y_hat = xXconv.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 
