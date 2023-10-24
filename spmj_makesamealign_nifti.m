function spmj_makesamealign_nifti(P,Q)
% function spmj_makesamealign_nifti(P,Q)
% This function uses the affine transformation matrix from the source image 
% and applies it to a 3d- or 4d- nifti file to give it the same affine transformation 
% matrix. It does not change the voxel data of those images. 
% It also deletes and <image>.mat files, where spm may have shuffled away
% different matrices for different slice of a 4d- nifti. 
% 
% Input:
%   P: name, vol-structure, or affine transformation matrix of the source image 
%   Q: Cell array of imagenames of the 4d-nifti files to apply the affine
%   matrix to. 
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
 