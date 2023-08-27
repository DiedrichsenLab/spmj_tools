function [beta,ResMS,y_hat,y_res,r2]=spmj_est_roi(Y,xX);
% function [beta,ResMS,y_hat,y_res,r2]=spmj_est_roi(Y,xX);
% INPUT:
%       Y:      N x P matrix of raw data 
%       xX:     Design structure for FIR estimation 
%       xXconv: Design structure for convolution (both from spmj_svdhrf_prep)
% OUTPUT:
%       beta:   Estimated beta (QxP) 
%       ResMS:  residual mean of sum of square (1xP)  
%       HRF:    the HRF that was used  (numFIRx1)
%       res:    full NxP residuals
% Joern Diedrichsen

KWY   = spm_filter(xX.K,xX.W*Y);       %- Filter the data
if (~isfield(xX,'pKX')) 
   xX.xKXs   = spm_sp('Set',spm_filter(xX.K,xX.W*xX.X));       % KWX
   xX.xKXs.X = full(xX.xKXs.X);
   xX.pKX    = spm_sp('x-',xX.xKXs);                        % projector
end; 
beta  = xX.pKX*KWY;                    %- Parameter estimates
P     = size(beta,2);
y_res   = spm_sp('r',xX.xKXs,KWY);       %-Residuals
df    = size(xX.X,1)-size(xX.X,2);     % Df 
ResSS = sum(y_res.^2); 
ResMS = ResSS./df;               %-Residual MS
r2    = 1-sum(KWY.^2)./ResSS; 

reg_interest=[xX.iH xX.iC]; 
y_hat = xX.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 

