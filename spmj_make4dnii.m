function spmj_make_4dnii(outname,P);
% function spmj_make_4dnii(outname,P);
% Compresses a series of 3Dnii files into a 4d-nifti
% INPUT: 
%   outname: file name of the final file 
%   P: Character array of the input files
% OUTPUT: 
%   is saved under output name 
if (nargin<2 || isempty(P))
    P=spm_select(inf,'image','to Transform');
end; 
V=spm_vol(P); 
for i=1:length(V)
    X=spm_read_vols(V(i)); 
    V(i).fname=(outname); 
    V(i).mat=V(1).mat;
    V(i).n=[i 1];
    spm_write_vol(V(i),X); 
    fprintf('.'); 
end;
fprintf('\n');