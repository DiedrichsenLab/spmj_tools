function [SPM,Yhat,Yres, p_opt]=spmj_fit_hrfparams_example()
% Example script for the adjustment of the hemodynamic response function
%  The standard parameters are: 
%        p(1) - delay of response (relative to onset)          6
%        p(2) - delay of undershoot (relative to onset)       16
%        p(3) - dispersion of response                         1
%        p(4) - dispersion of undershoot                       1
%        p(5) - ratio of response to undershoot                6
%        p(6) - onset {seconds}                                0
%        p(7) - length of kernel {seconds}                    32
my_cifti=cifti_read('y_raw.dtseries.nii');
brain_model = my_cifti.diminfo{1}.models{1};  % Get the first brain model 
% Yraw = double(my_cifti.cdata(brain_model.start:brain_model.start+brain_model.count-1,:)'); % Pick out data as a T x P matrix 
Yraw = double(my_cifti.cdata(:,:)'); % Pick out data as a T x P matrix 
clear my_cifti; % preserve memory 

load SPM; % loading SPM structure 

% Get the hemodynamic response in micro-time resolution
SPM.xBF.UNITS    = 'secs'; % units of the hrf
SPM.xBF.T        = 16; % microtime resolution: number of time samples per scan
SPM.xBF.dt       = SPM.xY.RT/SPM.xBF.T;  % Delta-t per sample 
SPM.xBF.T0       = 1; % first time bin
SPM.xBF.name     = 'fitted_hrf';
SPM.xBF.order    = 1;
SPM.xBF.Volterra = 1;  % volterra expansion order?
SPM.xBF.bf = spm_hrf(SPM.xBF.dt,[6 16 1 1 6 0 32]);
SPM.xBF.length = size(SPM.xBF.bf,1)*SPM.xBF.dt; % support in seconds 

% Reconvolve the design matrix with new HRF
SPM = spmj_glm_convolve(SPM);

% Restimate the betas and get predicted and residual response
[beta, Yhat, Yres] = spmj_glm_fit(SPM,Yraw);

% Diagnostic plots 
subplot(2,2,1);
t=[SPM.xBF.dt*SPM.xBF.T0:SPM.xBF.dt:SPM.xBF.length]; 
plot(t,SPM.xBF.bf);
drawline([0],'dir','horz','linestyle','-');
title('Basis Function(s)'); 

% First run, convolved design matrix, sum of regressors of interest
subplot(2,2,2);
X=SPM.xX.X(SPM.Sess(1).row,SPM.Sess(1).col);
plot(sum(X,2));
title('Run 1 - overall response');

subplot(2,2,3);
Yhat_run1 = mean(Yhat(SPM.Sess(1).row,:),2);
Yres_run1 = mean(Yres(SPM.Sess(1).row,:),2);
Yadj_run1 = Yres_run1+Yres_run1; 
t= SPM.Sess(1).row; 
plot(t,Yhat_run1,'b',t,Yadj_run1,'b:'); 
title('Run 1 - overall response');

% Get onset structure, cut-out the trials of choice, and plot evoked
% response
subplot(2,2,4);
D = spmj_get_ons_struct(SPM);
Yadj = Yres+Yres; 
pre = 10;
post = 20;
for i=1:size(D.block,1);
    D.y_adj(i,:)=cut(mean(Yadj,2),pre,round(D.ons(i))-1,post,'padding','nan')';
    D.y_hat(i,:)=cut(mean(Yhat,2),pre,round(D.ons(i))-1,post,'padding','nan')';
    D.y_res(i,:)=cut(mean(Yres,2),pre,round(D.ons(i))-1,post,'padding','nan')';
end;
                
T = getrow(D,mod(D.num,2)==1); % Get the first onset for each double 
traceplot([-pre:post],T.y_adj,'errorfcn','stderr'); % ,
hold on;
traceplot([-pre:post],T.y_hat,'linestyle','--',...
        'linewidth',3); % ,
drawline([-8 8 16],'dir','vert','linestyle',':');
drawline([0],'dir','vert','linestyle','--');
drawline([0],'dir','horz','linestyle','-');
hold off;
xlabel('TR');
ylabel('activation');


