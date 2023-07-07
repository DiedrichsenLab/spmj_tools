function [P,SPM,Yhat,Yres]=spmj_fit_hrf(SPM,Yraw,varargin);
% Fits the two-gamma hrf fucntion to extracted time 
% function [hrf,P,xBF]=spmj_fit_hrf(SPM,Y_raw,options);
% INPUT: 
%   SPM:    SPM-structure with model to be fitted 
%   Yraw:   raw time series for hrf fitting 
% OUTPUT: 
%   P:      parameters fitted 
%   SPM:    changed SPM
%   Yhat:   predicted time series 
%   Yres:   residual time series 
% Parameters that can be fitted: 
%   p(1) - delay of response (relative to onset)         6
%   p(2) - delay of undershoot (relative to onset)      16
%   p(3) - dispersion of response                        1
%   p(4) - dispersion of undershoot                      1
%   p(5) - ratio of response to undershoot               6
%   p(6) - onset (seconds)                               0
%   p(7) - duration of the event ?? 
% --------------------------------------------------------------
P0     = [6  16  1   1    6  0  10]';               % Default parameters for the SPM hrf
LB     = [1  9   0.2 0.2  3  -2 0.2]';     
UB     = [11 25  4   4    12 6  15]'; 
UB     = [20 35  6   6    12 6  15]'; 

fit    = [1:5 7]';                  % Which of these parameters should be fit, last on is duration of event 
prior  = [0.4 0.3 1 1 0.5 1 0.2]'/800; 

vararginoptions(varargin,{'P0','fit','prior'}); 

% Fill in the default duration of the event, if not given 
if (length(P0)<7)
     P0(7)=SPM.Sess(1).U(1).dur(1); 
end; 

T=size(Yraw,1);
nscans = SPM.nscan;

SPM.xBF.T        = 16;  % time bins per scan
SPM.xBF.T0       = 1;  % first time bin
SPM.xBF.UNITS    = 'scans'; % units of the 
SPM.xBF.name     = 'fitted_hrf';
SPM.xBF.length   = 32.125;  % support in seconds
SPM.xBF.order    = 1;
SPM.xBF.Volterra = 1;  % volterra expansion order?
SPM.xBF.dt       = SPM.xY.RT/16;

% Filter and prepare the data 
Y = spm_filter(SPM.xX.K,SPM.xX.W*Yraw);

% Reduce the parameters to the fitted values and expand again later 
th0=P0(fit); 
th=fmincon(@cost,th0,[],[],[],[],LB(fit),UB(fit),[],[],Y,SPM,fit,P0,prior(fit)); 
% th=fminsearch(@cost,th0,[],Y,SPM,fit,P0); 
P=P0;
P(fit)=th;

SPM.xBF.bf=spmj_hrf(SPM.xBF.dt,P(1:7));
SPM=spmj_fMRI_design_changeBF(SPM); 

% return predicted timeseries and residuals 
beta  = SPM.xX.pKX*Y;                    %-Parameter estimates
Yres  = spm_sp('r',SPM.xX.xKXs,Y);                        % get the 
reg_interest=[SPM.xX.iH SPM.xX.iC]; 
Yhat   = SPM.xX.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 

function err=cost(th,Y,SPM,fit,P0,prior);               % Returns accumulated fitting error 
P          = P0; 
P(fit)     = th; 
SPM.xBF.bf = spmj_hrf(SPM.xBF.dt,P(1:7));
SPM        = spmj_fMRI_design_changeBF(SPM); 
res        = spm_sp('r',SPM.xX.xKXs,Y);                        % get the 
err        = sum(sum(res.^2))/numel(res)+sum(prior.*((P(fit)-P0(fit)).^2));                       % Overall squared error 


        