function spmj_resample(PI,PO,factor)
% function spmj_resample(PI,PO,factor)
% Resamples a image with a new resolution 
% scales the mat-file and the dimesions of the image accordingly 
% For half resolution sample at factor 0.5
VI          = spm_vol(PI);
VO          = VI;
VO.fname    = deblank(PO);

% Calculate new dimensions and mat-matrix
if (length(factor)==1)
    factor=[factor factor factor];
end;
VO.dim(1:3)=round(VI.dim(1:3).*factor);
m=spm_matrix([0 0 0 0 0 0 1./factor]); 
VO.mat      = VI.mat*m;
hld=1;

VO = spm_create_vol(VO);
for x3 = 1:VO.dim(3),
        M  = inv(spm_matrix([0 0 -x3 0 0 0 1 1 1])*inv(VO.mat)*VI.mat);
        v  = spm_slice_vol(VI,M,VO.dim(1:2),hld);
        VO = spm_write_plane(VO,v,x3);
end;
spm_close_vol(VO);