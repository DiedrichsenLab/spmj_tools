function sd=smnooth_kernel3D_fast(d,s) 
% 3 dimensional convolution of an matrix
% FORMAT spmj_smooth(d, s)
% d  - data to be smoothed
% s  - [sx sy sz] Guassian filter width
%
% Tobias Wiestler 7/09

if (s(1)==0 & s(2)==0 & s(3)==0)
    sd= d; 
else
    s  = max(s,ones(size(s)));			% lower bound on FWHM
    s  = s/sqrt(8*log(2));				% FWHM -> Gaussian parameter

    x  = round(6*s(1)); x = [-x:x];
    y  = round(6*s(2)); y = [-y:y];
    z  = round(6*s(3)); z = [-z:z];
    x  = exp(-(x).^2/(2*(s(1)).^2));
    y  = exp(-(y).^2/(2*(s(2)).^2));
    z  = exp(-(z).^2/(2*(s(3)).^2));
    x  = x/sum(x);
    y  = y/sum(y);
    z  = z/sum(z);

    i  = (length(x) - 1)/2;
    j  = (length(y) - 1)/2;
    k  = (length(z) - 1)/2;

    sd= zeros(size(d));
    spm_conv_vol(d,sd,x,y,z,-[i,j,k]);
end