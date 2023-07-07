function spmj_make_t_contrast(SPM)
% function spmj_make_t_contrast(SPM)
xCon=SPM.xCon;
co=length(xCon);
j=1;
for i=SPM.xX.iC
    c=zeros(size(SPM.xX.X,2),1);
    c(i)=1;
    xCon(co+j)=spm_FcUtil('Set',SPM.xX.name{i}, 'T', 'c',c,SPM.xX.xKXs);
    j=j+1;
end;
SPM.xCon=xCon;
keyboard;
spm_contrasts(SPM,[co+1:co+j-1]);
keyboard;

    