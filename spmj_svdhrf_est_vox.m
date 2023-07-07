function [beta,ResSS,HRF,res,r2]=spmj_svdhrf_est_vox(Y,xX,xXconv);
% function [beta,ResSS,HRF,res,r2]=spmj_svdhrf_est_vox(Y,xX,xXconv);
% INPUT:
%       Y:      N x P matrix of raw data
%       xX:     Design structure for FIR estimation
%       xXconv: Design structure for convolution (both from spmj_svdhrf_prep)
% OUTPUT:
%       beta:   Estimated beta (QxP)
%       ResSS:  residual sum of square (1xP)
%       HRF:    the HRF that was used  (numFIRx1)
%       res:    full NxP residuals
%       r2:     R2 value for each voxel
% Joern Diedrichsen

KWY   = spm_filter(xX.K,xX.W*Y);       %- Filter the data
beta1  = xX.pKX*KWY;                    %- Parameter estimates
numbasis = max(xX.iCb);
P     = size(beta1,2);
beta  = zeros(size(xXconv.X,2),P); 
y_res   = zeros(size(xXconv.X,1),P); 
reg_interest=[xXconv.iH xXconv.iC]; 

for p=1:P
    B=reshape(beta1(xX.iCb>0,p),numbasis,sum(xX.iCb==1));
    [v,a]=svds(B,1);
    HRF(:,p)=v.*sign(mean(v(3:5)));
    
    % Now estimate the betas based on the estimated HRF
    i1              =  xXconv.iCc>0;  % get the columns to convolve
    a               =  conv2(xXconv.X(:,i1),HRF(:,p));
    xXconv.Xv(:,i1 )= a(1:size(xXconv.X,1),:);
    xXconv.xKXs     = spm_sp('Set',xXconv.Xv);
    xXconv.pKX      = spm_sp('x-',xXconv.xKXs);        % Pseudo inverse
    beta(:,p)       = xXconv.pKX*KWY(:,p);
    y_res(:,p)      = spm_sp('r',xXconv.xKXs,KWY(:,p));       %-Residuals
    y_hat(:,p)      = xXconv.xKXs.X(:,reg_interest)*beta(reg_interest,p); %- predicted values 
end;
ResSS          = sum(y_res.^2);                   %-Residual SS
r2             = 1-sum(KWY.^2)./ResSS;

