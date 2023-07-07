function err=spmj_hrflmw_costf(P,t,y,hrf)
% function y=spmj_hrflmw(t,P,hrf)
% Hemodynamic response function with parameters PS
% convolved with a boxcar function of latency, magnitude and width (1 s)
err=y-spmj_hrflmw(t,P,hrf);
