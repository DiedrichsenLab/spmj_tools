function SPM=spmj_setContrasts(C,names,nums,SPM,optimize); 
% function SPM=spmj_setContrasts(C,names,SPM,[optimize==1]); 
% Sets t-contrasts num to names and the the rows in C; 
% Key idea is to optimize the contrast weights to take into account the
% estimation variance. 
% This should make the result of a t-test using the contrast equivalent in
% the numerator to running a GLM with a collapsed design matrix X*C. 
% To achieve this, however, the intercept and other regressors need to be
% pulled through correctly. 
% INPUTS: 
%       C:       QxN contrast vector with Q being the number of contrasts and N
%                the number of task-related regressors 
%       names:   a Q-long cell array with names for the contrasts 
%       nums:    a Q-long vector that tells the routine wich contrasts in the
%                SPM structure need to be set 
%       SPM:     The SPM structure 
%       optimize:Set 1 to optimize, 0 for normal contrast (just use C)  
C=C(:,SPM.xX.iC);                           % Only get the columns that are task related. 
CA = blockdiag(C,eye(length(SPM.xX.iB)));  % Pull through the intercepts 
Z  = pinv(CA);                              % Reconstruct the design matrix that would underly that contrast 

X  = SPM.xX.xKXs.X; 
iV = X'*X;   % Inverse variance-covariance of the beta-weights from the first level 
if (nargin<5 || optimize==1) 
%     CA = (Z'*iV*Z)\(Z'*iV); 
    CA = pinv(Z'*iV*Z)*(Z'*iV); 
end; 

% Now set the contrast in the SPM structure 
for i=1:size(C,1)
%     SPM.xCon(nums)=spm_FcUtil('Set',names{i}, 'T', 'c',CA(i,:)',SPM.xX.xKXs);
    SPM.xCon(nums(i))=spm_FcUtil('Set',names{i}, 'T', 'c',CA(i,:)',SPM.xX.xKXs);
end; 
