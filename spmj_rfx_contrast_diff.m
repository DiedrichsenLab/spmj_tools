function SPM=spmj_rfx_contrast_diff(SPM);
% Makes all possible T-contrast of a basic model
% uses both the H and C matrices for contrast
% H-indexes effect of interest, C-confounds, B-block level
% Use this function to generate contrasts for between-group (2nd level)
% GLMs.

Contrasts=[SPM.xX.iH SPM.xX.iC SPM.xX.iB];
for i=1:length(SPM.xX.iH)
    c=zeros(size(Contrasts,2),1);
    c(Contrasts(i))=1;
    xCon(i)=spm_FcUtil('Set',SPM.xX.name{i},'T','c',c,SPM.xX.xKXs);
end;
% Now do the difference
if (length(SPM.xX.iH>1)) % More than 1 group: add difference contrast
    c=[-1 1]';
    xCon(i+1)=spm_FcUtil('Set',SPM.xX.name{i},'T','c',c,SPM.xX.xKXs);
end;
SPM.xCon=xCon;
SPM=spm_contrasts(SPM,1:length(xCon));
save SPM SPM;
