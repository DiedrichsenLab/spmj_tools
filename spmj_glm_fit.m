function [beta, Yhat, Yres]= spmj_glm_fit(SPM,Yraw)
% function beta, Yhat, Yres = spmj_glm_fit(SPM,Yraw)
% Estimates a new GLM on some raw time-series data. 
% Applies temporal filtering and whitening to the data before fitting. 
% INPUT: 
%   SPM: SPM-structure with convolved design matrix, filtering and
%        whitening matrices 
%   Yraw: TxP matrix of time series data. Usually from an ROI 
Y = spm_filter(SPM.xX.K,SPM.xX.W*Yraw);
beta  = SPM.xX.pKX*Y; %-Parameter estimates
Yres  = spm_sp('r',SPM.xX.xKXs,Y); % get the 
reg_interest=[SPM.xX.iH SPM.xX.iC]; 
Yhat   = SPM.xX.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 
