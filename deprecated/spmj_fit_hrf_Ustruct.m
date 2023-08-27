function P=spmj_fit_hrf_Ustruct(SPM,y,varargin);
% Fits the two-gamm function hrf to an adjusted time series of a 
% region 
% function [hrf,P,xBF]=spmj_fit_hrf(SPM,y_adj,options);
% xBF: basis function Structure
% --------------------------------------------------------------
model=1;

pre=3; 
post=8; 
fig=0;

vararginoptions(varargin,{'pre','post','fig'}); 

nscans = SPM.nscan;
t=[-pre:post]*SPM.xY.RT; 


T=size(y,1);
xBF.T=16;  % time bins per scan
xBF.T0=1;  % first time bin
xBF.UNITS='scans'; % units of the 
xBF.name='fitted_hrf';
xBF.length=32.125;  % support in seconds
xBF.order=1;
xBF.Volterra=1;  % volterra expansion order?
xBF.dt=SPM.xY.RT/16;

% Extract onsets from SPM
UU.ons=[];
UU.dur=[];
UU.u=[];

    offset=0;
for b=1:length(SPM.nscan)
    D=size(SPM.Sess(b).U(1).u,1);
    U(b).ons=[];U(b).u=spalloc(D,1,D);U(b).dur=[];
    Ub=spm_get_ons(SPM,b);
    for u=1:length(Ub);
        U(b).ons=[U(b).ons;Ub(u).ons];
        U(b).u=[U(b).u + Ub(u).u];
        U(b).dur=[U(b).dur;Ub(u).dur];
        U(b).name={'event'};
    end;
end;
    

P0=[6 16 1 1 6 0]';

P=fminsearch(@(p) cost(p,U,y),P0(1:model)); 
y_hat=predicted(P,U,y);

[D,session_onset] = spmj_get_ons_struct(SPM);
for i=1:length(D.block) 
    D.y_adj(i,:)=cut(y,pre,round(D.ons(i)),post,'padding','nan'); 
    D.y_hat(i,:)=cut(y_hat,pre,round(D.ons(i)),post,'padding','nan'); 
end; 

if (fig==1)
    traceplot(t,D.y_adj,'linestyle','-'); hold on;
    traceplot(t,D.y_hat,'linestyle',':');hold off;
end;


function e=cost(p,U,y);
y_hat=predicted(p,U,y);
r=y-y_hat;
e=r'*r;

function y_hat=predicted(P,U,y);
N=length(y);
bf=spm_hrf(2.015/16,P);
Xx=[];
for b=1:length(U) 
    X = spm_Volterra(U(b),bf,1);
    X = X([0:(146 - 1)]*16 + 1 + 32,:);
    Xx = blkdiag(Xx,X);
end;
X=[Xx ones(length(y),1)];
y_hat=X*inv(X'*X)*X'*y;
