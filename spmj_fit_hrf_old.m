function [hrf,P,xBF]=spmj_fit_hrf(RT,y,fig,options);
% Fits the two-gamm function hrf to actual (or fittted) responses 
% function [hrf,P,xBF]=spmj_fit_hrf(RT,y,fig);
% RT   - scan repeat time
% hrf  - hemodynamic response function
% p    - parameters of the response function (two gamma functions)
%							defaults
%							(seconds)
%	p(1) - delay of response (relative to onset)	   6
%	p(2) - delay of undershoot (relative to onset)    16
%	p(3) - dispersion of response			   1
%	p(4) - dispersion of undershoot			   1
%	p(5) - ratio of response to undershoot		   6
%	p(6) - onset (seconds)				   0
%	p(7) - length of kernel (seconds)		  32
%
% xBF: basis function Structure
% --------------------------------------------------------------
if (nargin<4)
    options='full';
end;

T=size(y,1);
le=RT*(T-1);
yscale=1;
t=([1:T]-1)*RT;
xBF.T=16;  % time bins per scan
xBF.T0=1;  % first time bin
xBF.UNITS='scans'; % units of the 
xBF.name='fitted_hrf';
xBF.length=32.125;  % support in seconds
xBF.order=1;
xBF.Volterra=1;  % volterra expansion order?
xBF.dt=1/8;


OPTIONS=optimset('lsqnonlin');
OPTIONS.MaxFunEvals=3000;
OPTIONS.TolFun=0.0000000001;
switch (options)
    case 'full'
        P0=[6 16 1 1 6 0 yscale];
        LB=[0 0  0 0 0 -inf 0];  
        [P]=lsqnonlin(@spmj_fit_hrf_costf,P0,LB,[],OPTIONS,RT,le,y);
        yscale=P(7);
        P(7)=le;
        hrf=spm_hrf(RT,P)*yscale;
        P(7)=xBF.length;
        xBF.bf=spmj_hrf(xBF.dt,P);
    case 'restricted'
        P0=[6 16 1 1 6 0 le];
        yscale=lsqnonlin(@spmj_fit_hrf_costf_res,yscale,0,[],OPTIONS,RT,y,P0);
        P=P0;
        hrf=spm_hrf(RT,P)*yscale;
        P(7)=xBF.length;
        xBF.bf=spmj_hrf(xBF.dt,P);
    case 'lmw' % Latency magnitude width
        LB=[-2 0 0.125];
        UB=[4 inf 30];
        
        P0={[0 20 0.2],[2 1 5],[-1 2 3]};
        for i=1:length(P0)
            [P{i},resnorm(i)]=lsqnonlin(@spmj_hrflmw_costf,P0{i},LB,UB,OPTIONS,t',y,spm_hrf(1/8));
        end;
        resnorm
        [dummy,indx]=min(resnorm);
        P=P{indx};
        hrf=spmj_hrflmw(t,P,spm_hrf(1/8));
end;    
if (nargin>2 & fig==1) 
    plot(t,y,'r',t,hrf,'r:');
end;
