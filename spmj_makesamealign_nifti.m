function M = spmj_makesamealign_nifti(P,Q)
% function M = spmj_makesamealign(P,Q)
% makes all images 
% have the same alignement 
% Does it together for a whole 4D image and deletes the 
% mat file! Use spmj_makesamealign if you want to modify only parts 
% of a NIFTI file 

if (nargin<1 || isempty(P)) 
    P=spm_select(1,'image','Select Image for orientation info'); 
end; 
if (nargin<2 || isempty(Q)) 
    Q=spm_select(inf,'image','Select Images to make equal'); 
end; 
if (isnumeric(P))
    VP.mat=P;
elseif (isstruct(P))
    VP=P; 
else     
    VP=spm_vol(P); 
end; 

fprintf('First Image:\n');
VP.mat 
k=0;
for i=1:size(Q,1)
    [dir,name,ext,num]=spm_fileparts(Q(i,:)); 
    N=nifti(fullfile(dir,[name ext])); 
    N.mat=VP.mat; 
    N.mat0=VP.mat; 
    create(N); 
    matfile=fullfile(dir,[name '.mat']); 
    if (exist(matfile,'file')) 
        delete(matfile); 
    end; 
end; 
fprintf('\n'); 
 