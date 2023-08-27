function [Z,v]=smooth_kernel3D(y,sigma);
% function y=smooth_kernel(y,sigma);
% runs a gaussian smoothing kernel over the data
% Does not intorduce shift of data 
if (length(sigma)==1) 
    sigma=[sigma sigma sigma];
end;
N=round(sigma*2)*2;
[X,Y,Z]=meshgrid([0:N(1)-1],[0:N(2)-1],[0:N(3)-1]);
D=sqrt(((X-N(1)/2)/sigma(1)).^2+((Y-N(2)/2)/sigma(2)).^2+((Z-N(3)/2)/sigma(3)).^2);
v=normpdf(D,0,1);v=v./sum(v(:));
Z=convn(y,v,'same');
