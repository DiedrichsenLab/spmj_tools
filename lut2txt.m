function X=lut2txt(fname)
% Transforms a color lut file to a text file 
f=fopen(fname); 
X=fread(f); 
numrows=length(X)/3; 
X=reshape(X,numrows,3); 
