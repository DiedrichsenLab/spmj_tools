function [I]=spmj_slice_vol(V,x,y,z);
% Takes a slice through one image 
% in X, Y or Z direction 
if (length(x)==1) 
    [Y,Z]=meshgrid(y,z);
    X=ones(size(Y))*x;Z=flipud(Z);
elseif (length(y)==1) 
    [Y,X]=meshgrid(x,z);Y=flipud(Y);
    Z=ones(size(X))*y;
elseif (length(z)==1) 
    [X,Y]=meshgrid(x,y);Y=flipud(Y);
    Z=ones(size(X))*z;
else 
    error ('one dimension must be singluar');
end;
[row,col]=size(X);
XYZ=[X(:) Y(:) Z(:) ones(row*col,1)]';

XYZn=inv(V.mat)*XYZ;
X=reshape(XYZn(1,:),row,col);
Y=reshape(XYZn(2,:),row,col);
Z=reshape(XYZn(3,:),row,col);
I=spm_sample_vol(V,X,Y,Z,0);

