function error=spmj_fit_hrf_costf(yscale,dt,y,P0);
yp=spm_hrf(dt,P0).*yscale;
error=y-yp;
