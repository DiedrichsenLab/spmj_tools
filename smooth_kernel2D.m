function [Z,v]=smooth_kernel2D(y,sigma);
% function y=smooth_kernel(y,sigma);
% runs a gaussian smoothing kernel over the data
% Does not intorduce shift of data 
N=round(sigma*5)*2;
[X,Y]=meshgrid([0:N-1],[0:N-1]);
v=normpdf(sqrt((X-(N)/2).^2+(Y-(N)/2).^2),0,sigma);
v=v./sum(sum(v));
Z=conv2(y,v,'same');
