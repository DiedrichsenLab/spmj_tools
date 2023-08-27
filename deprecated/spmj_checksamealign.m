function M = spmj_checksamealign(P,Q)
% function M = spmj_checksamealign(P,Q)
% check whether the images
% are all the same alignement 

if (nargin<1 || isempty(P)) 
    P=spm_select(1,'image','Select first Image'); 
end; 
if (nargin<2 || isempty(Q)) 
    Q=spm_select(inf,'image','Select Images to compare them to'); 
end; 
VP=spm_vol(P); 
VQ=spm_vol(Q); 

fprintf('First Image:\n');
VP.mat 
k=0;
for i=1:length(VQ)
    if any(any(abs(VP.mat-VQ(i).mat)>0.0001))
        fprintf('Image %s has different alignment:\n',VQ(i).fname);
        VQ(i).mat
        k=k+1; 
    else 
        fprintf('Image %s has same alignment:\n',VQ(i).fname);        
    end; 
end; 
