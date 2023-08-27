function [SPM]=spmj_make_contrast(SPM,C,Name);
% function [SPM]=coord4_make_contrast(SPM,C,Name);
% Makes the t-contrast specified with rows of C and Name 

for i=1:size(C,1)
    xCon(i)=spm_FcUtil('Set',Name{i}, 'T', 'c',C(i,:)',SPM.xX.xKXs);
end;
SPM.xCon=xCon;
SPM=spm_contrasts(SPM,[1:size(C,1)]);
save SPM SPM 
