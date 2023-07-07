function y=spmj_conv(x,basis)
% function spmj_conv(x,basis,varargin)
% Convolution of a vector x with a basis function basis 
% x: A column vector or matrix 
% basis: a basis function. If it has multiple columns, x is convolved with
%   each of these columns
% Joern Diedrichsen 2008
[xr,xc]=size(x);
[br,bc]=size(basis); 
if (xr<xc) 
    error ('x must be a column vector, or a N*p matrix'); 
end; 
if (br==1) 
    warning('basis should be a column vector'); 
end; 
y=[];
for i=1:xc
    for j=1:bc
        a=conv(x(:,i),basis(:,j)); 
        y=[y a(1:xr)]; 
    end;
end; 