function SPM=spmj_rfx_contrast(SPM);
% Makes all possible T-contrast of a basic model 
% uses both the H and C matrices for contrast
Contrasts=[SPM.xX.iH SPM.xX.iC SPM.xX.iB];
for i=1:length(Contrasts)
    c=zeros(size(Contrasts,2),1);
    c(Contrasts(i))=1;
    xCon(i)=spm_FcUtil('Set',SPM.xX.name{i},'T','c',c,SPM.xX.xKXs);
end;
SPM.xCon=xCon;
SPM=spm_contrasts(SPM,1:length(Contrasts));
save SPM SPM;
