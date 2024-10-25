function issame = spmj_check_sform_vs_qform(fname)
% function issame = spmj_check_sform_vs_qform(fname)
% Checks is qform and sform are the same 
% INPUT: 
%       fname: Nifti file name 
% RETURNS: 
%       issame: True if both are approximately the same 
eps = 1e-4; % Tolerance 
V=spm_vol(fname); 
sform = V.private.mat;
qform = V.private.mat0;
if sum(sum(abs(sform-qform))) < eps
    fprintf('Sform and Qform are approximately the same\n'); 
    issame = 1;  
else 
    fprintf('Sform and Qform are different:\n'); 
    sform
    qform
    issame = 0; 
end
    
    