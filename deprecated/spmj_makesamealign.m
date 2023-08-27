function M = spmj_makesamealign(P,Q)
% function M = spmj_makesamealign(P,Q)
% makes all images 
% have the same alignement 

if (nargin<1 || isempty(P)) 
    P=spm_select(1,'image','Select Image for orientation info'); 
end; 
if (nargin<2 || isempty(Q)) 
    Q=spm_select(inf,'image','Select Images to make equal'); 
end; 
VP=spm_vol(P); 
VQ=spm_vol(Q); 

fprintf('First Image:\n');
VP.mat 
k=0;
for i=1:length(VQ)
    X=spm_read_vols(VQ(i)); 
    VQ(i).mat=VP.mat; 
    spm_write_vol(VQ(i),X); 
    clc; 
    fprintf('%d',i); 
end; 
fprintf('\n'); 
 