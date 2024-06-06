function [SPM,Yhat,Yres, p_opt]=spmj_fit_hrfparams(SPM,Yraw)
% Fits the two-gamma hrf fucntion to extracted time 
% function [hrf,P,xBF]=spmj_fit_hrf(SPM,Y_raw,options);
% INPUT: 
%   SPM:    SPM-structure with model to be fitted 
%   Yraw:   raw time series for hrf fitting 
% OUTPUT: 
%   SPM:    changed SPM
%   Yhat:   predicted time series 
%   Yres:   residual time series 
%   p_opt: parameters fitted
% p    - parameters of the response function (two Gamma functions)
%                                                           defaults
%                                                          {seconds}
%        p(1) - delay of response (relative to onset)          6
%        p(2) - delay of undershoot (relative to onset)       16
%        p(3) - dispersion of response                         1
%        p(4) - dispersion of undershoot                       1
%        p(5) - ratio of response to undershoot                6
%        p(6) - onset {seconds}                                0
%        p(7) - length of kernel {seconds}                    32
% --------------------------------------------------------------


SPM.xBF.T        = 16; % time bins per scan
SPM.xBF.T0       = 1; % first time bin
SPM.xBF.UNITS    = 'secs'; % units of the 
SPM.xBF.name     = 'fitted_hrf';
% TODO
SPM.xBF.length   = 50.0625; % support in seconds 
SPM.xBF.order    = 1;
SPM.xBF.Volterra = 1;  % volterra expansion order?
SPM.xBF.dt       = SPM.xY.RT/16;

% Filter and prepare the data 
Y = spm_filter(SPM.xX.K,SPM.xX.W*Yraw);
% optimization options
option.TolX = 1E-10;
option.Display = 'off';
option.MaxFunEvals = 500;

%% Optimizing four parameters, p(1), p(2), p(6), p(7)
% options = optimoptions('fmincon','Algorithm','interior-point'); % run interior-point algorithm
% LB = [2 5 0 16]; UB = [10 20 5 40];
% p_opt = fmincon(@(p) cost(p, Y,SPM), p0, [1 -1 0 0;0 1 0 -1],zeros(2,1), [], [], LB, UB,[],options);
% SPM.xBF.bf = spm_hrf(SPM.xY.RT/16, [p_opt(1:2)' 1 1 6 p_opt(3:4)'], 16);  % replace the basis function with optimal hrf
%% Optimizing two parameters, p(1), p(2)
p_def = spm_get_defaults('stats.fmri.hrf');
options = optimoptions('fmincon','Algorithm','interior-point'); % run interior-point algorithm
LB = [2 5]; UB = [10 20];
p_opt = fmincon(@(p) cost(p, Y,SPM), p_def, [1 -1],0, [], [], LB, UB,[],options);

SPM.xBF.bf = spm_hrf(SPM.xY.RT/16, [p_opt(1:2) p_def(3:end)], 16);  % replace the basis function with optimal hrf


SPM = fMRI_design_changeBF(SPM); 

% return predicted timeseries and residuals 
beta  = SPM.xX.pKX*Y; %-Parameter estimates
Yres  = spm_sp('r',SPM.xX.xKXs,Y); % get the 
reg_interest=[SPM.xX.iH SPM.xX.iC]; 
Yhat   = SPM.xX.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 

function err=cost(p, Y,SPM)
% cost function to be minimized for parameter fitting
p_hrf = spm_get_defaults('stats.fmri.hrf');
p_hrf(1:2) = p; %% four parameters of interest (p(1), p(2), p(6), p(7))
% p_hrf([1:2 6 7]) = p; %% four parameters of interest (p(1), p(2), p(6), p(7))
SPM.xBF.bf = spm_hrf(SPM.xY.RT/SPM.xBF.T, p_hrf, SPM.xBF.T);
SPM = fMRI_design_changeBF(SPM);
res = spm_sp('r',SPM.xX.xKXs,Y); % get the residual
% constraints
err = sum(sum(res.^2))/numel(res);



function SPM = fMRI_design_changeBF(SPM)
% Re-convolves the SPM structure with a new basis function
% If second input argument is given, also changes duration of the event 
Xx    = [];
Xb    = [];
iCs   = [];     % number of the run 
iCc   = [];     % number of the condition
iCb   = [];     % number of the basis function 
iN    = [];     % Index for regressor of non-interest 
dur = 1;  %% 1 second duration, edited by SKim


numbasis=size(SPM.xBF.bf,2); 
numscan=length(SPM.nscan); 

for s = 1:numscan
    
    % number of scans for this session
    %----------------------------------------------------------------------
    k = SPM.nscan(s);

    
    % Get inputs, neuronal causes or stimulus functions U
    %------------------------------------------------------------------
    U=SPM.Sess(s).U;
    numcond=length(U); 
   
    % Convolve stimulus functions with basis functions
    %------------------------------------------------------------------
    [X,Xn,Fc] = spm_Volterra(U,SPM.xBF.bf,SPM.xBF.Volterra);
    
    % Resample regressors at acquisition times (32 bin offset)
    %------------------------------------------------------------------
    X = X((0:(k - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);
    
    % and orthogonalise (within trial type)
    %------------------------------------------------------------------
    for i = 1:length(Fc)
        X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
    end
    
    
    % get user specified regressors
    %==================================================================
    C     = SPM.Sess(s).C.C;
    numreg = size(C,2);
    X      = [X spm_detrend(C)];
    
    
    %-Session structure array
    %----------------------------------------------------------------------
    SPM.Sess(s).row    = size(Xx,1) + (1:k);
    SPM.Sess(s).col    = size(Xx,2) + (1:size(X,2));
       
    % Confounds: Session effects
    %==================================================================
    B      = ones(k,1);
    
    % append into Xx and Xb
    %======================================================================
    Xx    = blkdiag(Xx,X);
    Xb    = blkdiag(Xb,B);
    
    iCs   = [iCs ones(1,size(X,2)+numreg)*s];
    iCc   = [iCc kron([1:numcond],ones(1,numbasis)) zeros(1,numreg)]; 
    iCb   = [iCb kron(ones(1,numcond),[1:numbasis]) zeros(1,numreg)];
    iN    = [iN zeros(1,numcond*numbasis) ones(1,numreg)];
end

% finished
%--------------------------------------------------------------------------
SPM.xX.X      = [Xx Xb];
if isfield(SPM.xX,'W')
    SPM.xX.xKXs   = spm_sp('Set',spm_filter(SPM.xX.K,SPM.xX.W*SPM.xX.X));       % KWX
else
    SPM.xX.xKXs   = spm_sp('Set',spm_filter(SPM.xX.K,SPM.xX.X));       % KWX
end
SPM.xX.xKXs.X = full(SPM.xX.xKXs.X);
SPM.xX.pKX    = spm_sp('x-',SPM.xX.xKXs);                        % projector

SPM.xX.iC     = 1:size(Xx,2);
SPM.xX.iB     = (1:size(Xb,2)) + size(Xx,2);
% Indices 
SPM.xX.iCs    = [iCs [1:numscan]];   % number of the session
SPM.xX.iCc    = [iCc zeros(1,numscan)];   % number of the condition
SPM.xX.iCb    = [iCb zeros(1,numscan)];  % number of the basis function
SPM.xX.iN     = [iN  ones(1,numscan)*2];   % Flag for the regressors of no interest (+ intercept)

        