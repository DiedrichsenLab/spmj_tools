function spmj_segmentation_VBQ(anatomical_image, varargin)
%Segmentation for VBQ images (N. Weiskopf has tweaked unified segmentation for use with R1/MT images
%to apply UNICORT correction for RF transmit inhomogeneities to R1 maps estimated from dual angle FLASH
%exps.  The correction is primarily based on the new segment in SPM12b.):
%Sheena Waters-Metenier May2013

%***If any of the SPM12b segs fail, check that the matlabbatch extension (J) hasn't been renamed in your update***

%anatomical_image:  specify image you want to use for the initial affine registration with tissue probability map
%varargin:  specify whatever additional image(s) you want to use for mutual information (i.e. multi-modal/multi-spectral segmentation)

%Example usage:
%spmj_segmentation(ana_file);
%spmj_segmentation(ana_file, 'segOpt', 2, 'multimodal', {MM1_file, MM2_file, MM3_file});
%e.g. ana_file = 's01_1_MT.nii' and MM1_file = s01_1_R1.nii, MM2_file = s01_1_R2s.nii, MM3_file = s01_1_PDw.nii

%--------------------------------------------------------------------------
%To debug:
%Even when attempting via the gui, case 1 and 2 segmentation types yields 
%error 'index exceeds matrix dimensions'.
%Case 3 (old seg with modifed VBQ settings) does work.
%--------------------------------------------------------------------------

segOpt=1; % Default segmentation option: VBQ segmentation for DARTEL
multimodal={};
vararginoptions(varargin,{'segOpt','multimodal'});

spmVer=spm('Ver');
switch (segOpt)
    case 1                                                                 %VBQ Segmentation (for use with DARTEL, i.e. no norm applied)
        if ~strcmp(spmVer,'SPM12b')
            error('Only available in SPM12!');
        end;
                
        J.struct.s_vols     =  {[anatomical_image,',1']};                  %Specify the structural vol for unified segmentation.
        J.struct.biasreg    =  1;                                          %Select:  0=none, e.light=0.00001, v.light=0.0001, light=0.001, medium=0.01, heavy=0.1, v.heavy=1. e.heavy=10
        J.struct.biasfwhm   =  60;                                         %Select:  30,40,50,60,70,80,90,100,110,120,1230,140,150,None
        J.struct.write      =  [0 0];                                      %Select:  0=Save nothing.  1=save bias corrected (first) or bias field (second)
        
        spmLocation= which('SPM');

        
        %You can do DARTEL-importation at this stage by saving
        %DARTEL-imported image (r) in addition to the native.
        for i=1:6
            J.tissue(i).tpm =    {fullfile(spmLocation(1:end-6), 'tpm', ['TPM.nii,' num2str(i)])};  %Tissue probability map
            if i == [1 || 2 || 3 || 6]
                J.tissue(i).ngaus = 2;                                     %Select: 1 through 8 or 9 nonparametric   -->   Model here assumes that the intensity distribution of each tissue class may not be Gaussian and assigns probabilties according to these non-Gaussian distributions.  Gaussians are typically 2 for GM, WM, and CSF and 3 for bone, four for other soft tissues, and two for air.
            elseif i==4
                J.tissue(i).ngaus = 3;
            elseif i==5
                J.tissue(i).ngaus = 4;
            end;
            if i == [4 || 5 || 6]
                J.tissue(i).native = [0 0];                                %Select:  [1 0]=native, [0 1]=DARTEL imported, [1 1]=Native+DARTEL Imported   -->    Native space option allows you to generate tissue class image (c*) that is in alignment with the original.  It can also be used for importing into a form that can be used with DARTEL (rc*)
            else
                J.tissue(i).native = [1 0];
            end;
            J.tissue(i).warped = [0 0];                                    %Select:  [0 0]=nothing, [1 0]=modulated, [0 1]=unmodulated, [1 1]=modulated + unmodulated   -->    Modulation is to compensate for the effects of spatial normalsation.  When warping a series of images to match a template, it is inevitable that vol diff will be introduced into the warped images (e.g. if subj's temporal lobe is half that of the template, then its vol will be doubled during spat norm, which will double the num of voxels labelled as GM).  In order to remove this confound, the spat norm GM is adjusted by multiplying by it's relative volume pre- and post-warping.  If warping results in a region doubling its vol, then the correction will halve the intensity of the tissue label.  This procedure, essentially, has the effect of preserving the total amount of GM signal int he normalised partitions.
        end;
        
        J.warp.mrf       =   0;                                            %Select:  0=no cleanup, 1=clean
        J.warp.reg       =   4;                                            %Select:  if you select more regularisation, more smoothness (4 is the default, but this might be too much??  usually 1 is recommended).
        J.warp.affreg    =   'mni';                                        %Select: 'no affine reg', 'ICBM Euro'='mni', 'ICBM Asian', 'avg-sized template', or 'no regularisation'.
        J.warp.samp      =   3;                                            %Smaller vals are more accu.
        J.warp.write     =   [0 0];                                        %Select:  [0 0]=none, [1 0]=inverse, [0 1]=forward, [1 1]=inverse and forward.
        
        matlabbatch{1}.spm.tools.VBQ.proc.preproc8NoNorm= J;               %This is likely to change often.
    case 2                                                                 %VBQ Segmentation with normalisation
        if ~strcmp(spmVer,'SPM12b')
            error('Only available in SPM12!');
        end;
        %strcmp(spmVer,'SPM12b') && segOpt==3                              
        %This is *VBQ* new segmentation
        %Similar to Unified Segmentation, except for:
        %  (i)   slightly different treatment of the mixing proportions
        % (ii)   the use of an improved registration model
        %(iii)   the ability to use multi-spectral data
        % (iv)   an extended set of tissue probability maps, which allows a different treatment of voxels outside the brain
        %Note that the default settings are different, probably optimised for MT data.
        %See more complete descriptions of all variables in segOption 2 above.
        %Raw data must be entered in this order MT, PD, T1, B1, B0.
        
        %******************************************************************
        %IMPORTANT usage note for multi-modal:  
        %Input the anatomical image as usual (even though it will not be 
        %used) but also use it as the first image in your set of multimodal 
        %images.  SPM will read the first image in the set as the one for 
        %affine alignment!!! 
        %******************************************************************
        
        J.subjc(1).output.outdir = cd;                                     %Session 1 or 2 individual subj folder

        J.subjc(1).struct.s_vols    =   {[anatomical_image,',1']};         %Specify the volumes:  MT or T1 images for unified segmentation
        J.subjc(1).struct.biasreg   =   0.0001;                            %Select:  0=none, e.light=0.00001, v.light=0.0001, light=0.001, medium=0.01, heavy=0.1, v.heavy=1. e.heavy=10
        J.subjc(1).struct.biasfwhm  =   60;                                %Select:  30,40,50,60,70,80,90,100,110,120,1230,140,150,None
        J.subjc(1).struct.write     =   [0 0];                             %Selet:  0=Save nothing.  1=save bias corrected (first) or bias field (second)
        
        for i=1:length(multimodal)
            J.subjc(1).maps.mp_vols = {[multimodal{i},',1']};              %Specify the parameter maps for processing (MT, R2*, FA, etc).  Select all other images than the ones you're using for unified segmentation???
        end;
        
        spmLocation= which('SPM');
        
        for i=1:6
            J.tissue(i).tpm =    {fullfile(spmLocation(1:end-6), 'tpm', ['TPM.nii,' num2str(i)])};  %Tissue probability map
            if i == [1 || 2 || 3 || 6]
                J.tissue(i).ngaus = 2;                                     %Select: 1 through 8 or 9 nonparametric   -->   Model here assumes that the intensity distribution of each tissue class may not be Gaussian and assigns probabilties according to these non-Gaussian distributions.  Gaussians are typically 2 for GM, WM, and CSF and 3 for bone, four for other soft tissues, and two for air.
            elseif i==4
                J.tissue(i).ngaus = 3;
            elseif i==5
                J.tissue(i).ngaus = 4;
            end;
            if i == [4 || 5 || 6]
                J.tissue(i).native = [0 0];                                %Select:  [1 0]=native, [0 1]=DARTEL imported, [1 1]=Native+DARTEL Imported   -->    Native space option allows you to generate tissue class image (c*) that is in alignment with the original.  It can also be used for importing into a form that can be used with DARTEL (rc*)
            else
                J.tissue(i).native = [1 0];
            end;
            J.tissue(i).warped = [0 0];                                    %Select:  [0 0]=nothing, [1 0]=modulated, [0 1]=unmodulated, [1 1]=modulated + unmodulated   -->    Modulation is to compensate for the effects of spatial normalsation.  When warping a series of images to match a template, it is inevitable that vol diff will be introduced into the warped images (e.g. if subj's temporal lobe is half that of the template, then its vol will be doubled during spat norm, which will double the num of voxels labelled as GM).  In order to remove this confound, the spat norm GM is adjusted by multiplying by it's relative volume pre- and post-warping.  If warping results in a region doubling its vol, then the correction will halve the intensity of the tissue label.  This procedure, essentially, has the effect of preserving the total amount of GM signal int he normalised partitions.
        end;
        
        J.warp.mrf     =   0;                                              %Select:  0=no cleanup, 1=clean
        J.warp.reg     =   4;                                              %More regularisation, more smoothness (4 is the default, but this might be too much??  usually 1 is recommended).
        J.warp.affreg  =   'mni';                                          %Select 'no affine reg', 'ICBM Euro'='mni', 'ICBM Asian', 'avg-sized template', or 'no regularisation'.
        J.warp.samp    =   3;                                              %Select:  Smaller vals are more accurate.
        J.warp.write   =   [0 1];                                          %Select:  [0 0]=none, [1 0]=inverse, [0 1]=forward, [1 1]=inverse and forward.
        J.fwhm         =   [6 6 6];                                        %Specify:  smoothing kernel (must do this for modulated data; otherwise, you get aliaising)
        
        matlabbatch{1}.spm.tools.VBQ.proc.preproc8= J;                     %This is likely to change often.
    case 3                                                                 %This uses the OLD SEGMENT of SPM12b but with settings a bit more optimised for VBQ data
        J.data = {[anatomical_image,',1']};
        J.output.GM =  [0 0 1];
        J.output.WM =  [0 0 1];
        J.output.CSF = [0 0 1];
        J.output.biascor = 1;                                              
        J.output.cleanup = 0;
        spmLocation= which('SPM');
        J.opts.ngaus = [2 2 2 4];
        J.opts.regtype = 'mni';
        J.opts.warpreg = 4;                                                %instead of 1
        J.opts.warpco = 25;
        J.opts.biasreg = 1;                                                %heavy bias reg
        J.opts.biasfwhm = 60;
        J.opts.samp = 3;
        J.opts.msk = {''};
        
        if (strcmp(spmVer,'SPM12b'))
            J.opts.tpm = {
                fullfile(spmLocation(1:end-6),'toolbox', 'OldSeg', 'grey.nii')
                fullfile(spmLocation(1:end-6),'toolbox', 'OldSeg', 'white.nii')
                fullfile(spmLocation(1:end-6),'toolbox', 'OldSeg', 'csf.nii')
                };
            matlabbatch{1}.spm.tools.oldseg= J;
        else %e.g. SPM8
            J.opts.tpm = {
                fullfile(spmLocation(1:end-6),'tpm', 'grey.nii')
                fullfile(spmLocation(1:end-6),'tpm', 'white.nii')
                fullfile(spmLocation(1:end-6),'tpm', 'csf.nii')
                };
            matlabbatch{1}.spm.spatial.preproc= J;
        end;
end;


spm_jobman('run',matlabbatch);
