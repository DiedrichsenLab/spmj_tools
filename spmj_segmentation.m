function spmj_segmentation(anatomical_image, varargin)
%Run old or new segmentation (SPM12) with this function or old segmentation with SPM8 (see spmj_segmentation_VBQ for VBQ version).
%Sheena Waters-Metenier May2013

%***If any of the SPM12b segs fail, check that the matlabbatch extension (J) hasn't been renamed in your update***

%anatomical_image:  specify image you want to use for the initial affine registration with tissue probability map
%varargin:  specify whatever additional image(s) you want to use for mutual information (i.e. multi-modal/multi-spectral segmentation)

%Example usage:
%spmj_segmentation(ana_file);
%spmj_segmentation(ana_file, 'segOpt', 2, 'multimodal', {MM1_file, MM2_file, MM3_file});
%e.g. ana_file = 's01_1_MT.nii' and MM1_file = s01_1_R1.nii, MM2_file = s01_1_R2s.nii, MM3_file = s01_1_PDw.nii


segOpt=1; %Default segmentation option (old segmentation)
multimodal={};
vararginoptions(varargin,{'segOpt', 'multimodal'});

spmVer=spm('Ver');
switch (segOpt)
    case 1                                                                 %OLD SEGMENTATION
        J.data = {[anatomical_image,',1']};
        J.output.GM =  [0 0 1];
        J.output.WM =  [0 0 1];
        J.output.CSF = [0 0 1];
        J.output.biascor = 1;
        J.output.cleanup = 0;
        spmLocation= which('SPM');
        J.opts.ngaus = [2 2 2 4];
        J.opts.regtype = 'mni';
        J.opts.warpreg = 1;
        J.opts.warpco = 25;
        J.opts.biasreg = 0.0001;
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
    case 2                                                                 %NEW SEGMENTATION (with an option for multispectral data)    %For new segmentation, c1-c3 are as below and c4=skull and c5=meninges (c6s not saved)
        if ~strcmp(spmVer,'SPM12b')
            error('SegOption 2 is only available in SPM12');
        end;
        
        %You can do multi-spectral classification if you have data with different contrasts (but only the first channel will be used for the initial affine registration with tissue probability maps-->this is why I call it anatomical_image and name the other images something different (multispec_image).
        %NOTE:  if you use multiple channels, subj must be specified in the same order and all data must be in register (same position, size, voxel dimensions).
        
        spmLocation= which('spm');
        
        %******************************************************************
        %IMPORTANT usage note for multi-modal:  
        %Input the anatomical image as usual (even though it will not be 
        %used) but also use it as the first image in your set of multimodal 
        %images.  SPM will read the first image in the set as the one for 
        %affine alignment!!! 
        %******************************************************************
        
        J.channel(1).vols       =   {[anatomical_image,',1']};
        J.channel(1).biasreg    =   0.001;
        J.channel(1).biasfwhm   =   60;
        J.channel(1).write      =   [0 0];
        
        for i= 1:length(multimodal) %(+1 doesn't work)
            J.channel(i).vols       =   {[multimodal{i},',1']};
            J.channel(i).biasreg    =   0.001;                             %Select:  0=none, e.light=0.00001, v.light=0.0001, light=0.001, medium=0.01, heavy=0.1, v.heavy=1. e.heavy=10      MR images are often have 'bias'--that is, they are corrupted by a smooth, spatially varying artifact that modulates image intensity.  Bias impedes automatical processing.  Soo.. SPM models bias.  If your data is virtually no intensity uniformity, you will want to penalise large values for intensity non-uniformity parameters.  If data is more uniform, do not model bias.
            J.channel(i).biasfwhm   =   60;                                %Select:  30,40,50,60,70,80,90,100,110,120,1230,140,150,None        If your intensity non-uniformity is smooth, then select large FWHM-->this prevents algorithm from modelling out intensity variation due to different tissue types.  Model is i.i.d. Gaussian noise that has been smoothed by X amt, before taking the exponential.
            J.channel(i).write      =   [0 0];                             %Select:  0=Save nothing.  1=save bias corrected (first) or bias field (second)
        end;
        
        for i=1:6
            J.tissue(i).tpm =    {fullfile(spmLocation(1:end-6), 'tpm', ['TPM.nii,' num2str(i)])};  %Tissue probability map
            if i==1 || 2
                J.tissue(i).ngaus = 1;                                     %Select: 1 through 8 or 9 nonparametric   -->   Model here assumes that the intensity distribution of each tissue class may not be Gaussian and assigns probabilties according to these non-Gaussian distributions.  Gaussians are typically 2 for GM, WM, and CSF and 3 for bone, four for other soft tissues, and two for air.
            elseif  i==3 || 6
                J.tissue(i).ngaus = 2;
            elseif i==4
                J.tissue(i).ngaus = 3;
            elseif i==5
                J.tissue(i).ngaus = 4;
            end;
            if i==6                                                        
                J.tissue(i).native = [0 0];                                %Select:  [1 0]=native, [0 1]=DARTEL imported, [1 1]=Native+DARTEL Imported   -->    Native space option allows you to generate tissue class image (c*) that is in alignment with the original.  It can also be used for importing into a form that can be used with DARTEL (rc*)
            else
                J.tissue(i).native = [1 0];
            end;
            J.tissue(i).warped = [0 0];                                    %Select:  [0 0]=nothing, [1 0]=modulated, [0 1]=unmodulated, [1 1]=modulated + unmodulated   -->    Modulation is to compensate for the effects of spatial normalsation.  When warping a series of images to match a template, it is inevitable that vol diff will be introduced into the warped images (e.g. if subj's temporal lobe is half that of the template, then its vol will be doubled during spat norm, which will double the num of voxels labelled as GM).  In order to remove this confound, the spat norm GM is adjusted by multiplying by it's relative volume pre- and post-warping.  If warping results in a region doubling its vol, then the correction will halve the intensity of the tissue label.  This procedure, essentially, has the effect of preserving the total amount of GM signal int he normalised partitions.
        end;
        
        %SPM change note 1:  warped data here are not scaled by Jacobian determinants when generated the modulated data.  Rather, the original voxels are projected into their new location in the warped images. This exactly preserved tissue count, but has the effect of introducing aliasing artifacts--esp if the original data are at at lower res than the warped images.
        %SPM change note 2:  unmodulated data are also dealt with slightly different.  The projected data are corrected using a kind of smoothing procedure.  This is not done exactly as it should be (computationally speaking).  It also (imperfectly) extrapolates the warped tissue class images beyond the range of the original data.
        
        J.warp.mrf     =   1;                                                  %Select:  0=no cleanup, 1=clean    -->   Set the strength of the simple Markov Random Field.
        J.warp.cleanup =   1;                                                  %Select:  o=no clean, 1=light clean, 2=thorough clean     -->      This uses a crude routine for extracting brain from segmented images.  It begins by taking WM and eroding it a couple of times to get rid of any odd voxels.  The alogorithm continues to do conditional dilations for several iterations, where the condition is based upon GM or WM being present.  This identified region is then used to clean up the GM and WM partitions.  Note that the fluid class will also be cleaned, such that aqueous and virteos humour in the eyeballs (etc--except not CSF of course) will be removed.    %Note:  if pieces of your brains are chopped out, disable cleanup.
        J.warp.reg     =   [0 0.001 0.5 0.05 0.2];                             %Use pre-selected   -->   Warping regularisation:  registration involves min 2 terms:  (1) similarity msr b/t images (e.g. mean sq diff) and (2) roughness.  Roughness is the sum of:  (1) absolute displacements, (2) membrane energy, (3) bending energy, (4/5) linear elasticity regularisation, and (6) divergence.  (1)Absolute displacements must be penalised by a tiny amt. (2) The membrane energy of the deformation is penalised by a small amount (this penalises the sum of sq of the derivatives of the velocity field (the sum of sq of Jacobian tensors), (3) the bending energy penalises the sum of sq of the 2nd derivatives of the velocity, (4) linera elasticity reg encompases 2 elements (4 and 5).  The first (mu) is similar to that of linear elasticity, except it penalises the sum of sq of Jacobian tensors after they have been made symmetric by averaged with the transpose.  The term penalises length changes, with penalising rotations. (6) this weight denotes how much to penalise changes to the divergence of the velocities (lambda).  The divergence is a parameter of the rate of volumetric expansion or contraction.      %Regularisation determines the tradeoff between the 6 terms above.  More reg-->smoother deformations (smoothness is determined by bending energy of the deformations.
        J.warp.affreg  =   'mni';                                              %Select: 'no affine reg', 'ICBM Euro'=mni, 'ICBM Asian', 'avg-sized template', or 'no regularisation'.     -->    Affine regularisation:  this is a local optimisation and it needs reasonable starting estimates.  Img should be aligned already.  A mutual information affine reg is conducted with the tissue prob maps (note this does not model intensity non-uniformity), meaning that if the procedure is to be intialised with the affine reg, then the data should not be (too) corrupted with this artifact (you need to manually reposition if too much artifact).
        J.warp.fwhm    =   [0 0 0];                                            %Normally 0 for MRI.  This accounts for correlations between neighboring voxels (smoother data have more correlations; hence cannot be used with this).
        J.warp.samp    =   3;                                                  %The approx distance b/t sample pointed when estimating model params.  Smaller val use more of the dat, but this is slower.  Compromise between processing speed and accuracy.
        J.warp.write   =   [0 0];                                              %Select:  [0 0]=none, [1 0]=inverse, [0 1]=forward, [1 1]=inverse and forward   -->  Use with deformations utility.  For spat norm imgs to MNI space, you will need the forward deformation.
        
        matlabbatch{1}.spm.spatial.preproc= J;    
end;


spm_jobman('run',matlabbatch);
