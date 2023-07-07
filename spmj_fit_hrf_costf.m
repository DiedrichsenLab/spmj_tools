function error=spmj_fit_hrf_costf(P,dt,le,y);
yscale=P(7);
P(7)=le;
yp=spm_hrf(dt,P).*yscale;
error=y-yp;
