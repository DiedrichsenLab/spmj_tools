function spmj_bias_correct(Scans,varargin)
% function spmj_bias_correct(Scans,varargin)
% Runs a bias correction on a set of functional images
% Using the new segmentation tool
% INPUT:
%   Scans:      cell array on scans to biascorrect
% VARARGIN:
%   'segm', (0/1): 0 for skipping segmentation (1 is default)
% OUTPUT:
%   all scans will be written with prefix b
% New version for comaptibility with SPM12 
segm=1;
spmVer=spm('Ver');

% Estimate the bias field using new segmentation toolbox
vararginoptions(varargin,'segm');
spm_Dir= fileparts(which('spm'));
if (segm)
    global defaults
    J.channel.vols = {Scans{1}};
    J.channel.biasreg = 0.0001;
    J.channel.biasfwhm = 60;
    J.channel.write = [1 0];
    J.tissue(1).ngaus = 2;
    J.tissue(1).native = [1 0];
    J.tissue(1).warped = [0 0];
    J.tissue(2).ngaus = 2;
    J.tissue(2).native = [1 0];
    J.tissue(2).warped = [0 0];
    J.tissue(3).ngaus = 2;
    J.tissue(3).native = [1 0];
    J.tissue(3).warped = [0 0];
    J.tissue(4).ngaus = 3;
    J.tissue(4).native = [1 0];
    J.tissue(4).warped = [0 0];
    J.tissue(5).ngaus = 4;
    J.tissue(5).native = [1 0];
    J.tissue(5).warped = [0 0];
    J.tissue(6).ngaus = 2;
    J.tissue(6).native = [0 0];
    J.tissue(6).warped = [0 0];
    J.warp.affreg = 'mni';
    J.warp.samp = 3;
    J.warp.write = [0 0];
    switch(spmVer)
        case 'SPM12'
            J.tissue(1).tpm = {[spm_Dir '/tpm/TPM.nii,1']};
            J.tissue(2).tpm = {[spm_Dir '/tpm/TPM.nii,2']};
            J.tissue(3).tpm = {[spm_Dir '/tpm/TPM.nii,3']};
            J.tissue(4).tpm = {[spm_Dir '/tpm/TPM.nii,4']};
            J.tissue(5).tpm = {[spm_Dir '/tpm/TPM.nii,5']};
            J.tissue(6).tpm = {[spm_Dir '/tpm/TPM.nii,6']};
            J.warp.mrf = 1;
            J.warp.cleanup = 1;
            J.warp.reg = [0 0.001 0.5 0.05 0.2];
            J.warp.fwhm = 0;
            matlabbatch{1}.spm.spatial.preproc=J;
        case 'SPM8'
            J.tissue(1).tpm = {[spm_Dir '/toolbox/Seg/TPM.nii,1']};
            J.tissue(2).tpm = {[spm_Dir '/toolbox/Seg/TPM.nii,2']};
            J.tissue(3).tpm = {[spm_Dir '/toolbox/Seg/TPM.nii,3']};
            J.tissue(4).tpm = {[spm_Dir '/toolbox/Seg/TPM.nii,4']};
            J.tissue(5).tpm = {[spm_Dir '/toolbox/Seg/TPM.nii,5']};
            J.tissue(6).tpm = {[spm_Dir '/toolbox/Seg/TPM.nii,6']};
            J.warp.reg = 4;
            matlabbatch{1}.spm.tools.preproc8=J;
    end;
    spm_jobman('run',matlabbatch);
end;

[dir,name,ext,num]=spm_fileparts(Scans{1});
P=fullfile(dir,['BiasField_' name ext]);
Vin(1)=spm_vol(P);
for s=1:length(Scans)
    Vin(2)= spm_vol(Scans{s});
    Vout = Vin(2);
    [dir,name,ext,num]=spm_fileparts(Vout.fname);
    Vout.fname=fullfile(dir,['b' name ext num]);
    spm_imcalc(Vin,Vout,'i1.*i2');
end;
