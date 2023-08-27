function spmj_make3dnii;
% Transforms all images into 3d-nii's 
P=spm_select(inf,'image','to Transform');
V=spm_vol(P); 
for i=1:length(V)
    X=spm_read_vols(V(i)); 
    [pth,nam,ext,num] = spm_fileparts(V(i).fname);
    V(i).fname=fullfile(pth,[nam '.nii']); 
    spm_write_vol(V(i),X); 
    fprintf('.'); 
end; 
fprintf('\n'); 
