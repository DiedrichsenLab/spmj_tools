function spmj_makefieldmap(fmap_dir,magnitude_img, phasediff_img, varargin)
% function spmj_makefieldmap(varargin)
%   fmap_dir: Directory that holds the fmap and will be used for output
%   magnitude_img: Name of magnitude image 
%   phasediff_img: name of phase difference image 
% VARARGINOPTIONS 
%   func_dir: Directory that holds the unziped imaging_data 
%   epi_files: List of EPI-run files 
%   et1: TE1 
%   et2: TE2
%   tert: Total time 
% Prefix of the EPI data:
prefix = '';

% Which image (volume) to unwrap. This is mainly for the quality control 
% of field map correction. The corrected image will be saved with
% the prefix 'u'. 
image = 1;

% Directory of the EPI data:
func_dir = '';

% Name of the files to align to 
epi_files = {};

% fieldmap parameters:
et1 = 4.92;
et2 = 7.38;
tert = 0.7 * 90 * 1/2;
phase_encode = -1; % AP: -1 PA: +1 

% Handling the input arguments:
vararginoptions(varargin,{'func_dir','epi_files','et1', 'et2', 'tert','phase_encode'});

% Directory of the spm toolbox:
spm_dir = fileparts(which('spm'));

% displaying whether using 3D files or not
% disp(['Using 3D: ' num2str(use3D)])


% Parameters for the SPM field map toolbox:

%_______DEFAULTS Values_________________________________
% Echo times for magnitude images. Can be found in the meta data of the images. Default parameter for CFMM 3T scanner:
J.defaults.defaultsval.et = [et1 et2];

% If masking brain is on, the magnitude image is used to mask the brain:
J.defaults.defaultsval.maskbrain = 1;

% Phase-enocding direction. Should be available in the imaging sequence specifications. 
% Spm defines Anterior to Posterior as -1 and Posterior to Anterior as +1:
J.defaults.defaultsval.blipdir = phase_encode;

% Total EPI readout time = echo spacing (in ms) * base resolution (also knows as number of echos).  
% If you use GRAPPA acceleration, you need to divide the total number of echos by two:
J.defaults.defaultsval.tert = tert;

% EPI fieldmap or non-EPI. In CFMM 3T gradient echo imaging is used for
% field map acquisition which is a non-EPI sequence:
J.defaults.defaultsval.epifm = 0;

% Whether to apply Jacobian Modulation:
J.defaults.defaultsval.ajm = 0;

% DEFAULTS Values - uflags: for phase unwrapping and field map processing
% The unwrapping method:
J.defaults.defaultsval.uflags.method = 'Mark3D';

% FWHM of the Gaussian filter used to implement weighted smoothing of
% unwrapped images:
J.defaults.defaultsval.uflags.fwhm = 10;

% Size of padding kernel if required:
J.defaults.defaultsval.uflags.pad = 0;

% Normal or weighted smoothing:
J.defaults.defaultsval.uflags.ws = 1;

% DEFAULTS Values - mflags: for segmentation and creation of the brain mask
% Template file for segmentation to create brain mask:
J.defaults.defaultsval.mflags.template = {fullfile(spm_dir,'canonical','avg152T1.nii')};

% FWHM of Gaussian filter for smoothing brain mask:
J.defaults.defaultsval.mflags.fwhm = 5;

% Number of erosions used to create brain mask:
J.defaults.defaultsval.mflags.nerode = 2;

% Number of dilations used to create brain mask:
J.defaults.defaultsval.mflags.ndilate = 4;

% Threshold used to create brain mask from segmented data:
J.defaults.defaultsval.mflags.thresh = 0.5;

% value used in the segmentation. A larger value helps the segmentation to 
% converge:
J.defaults.defaultsval.mflags.reg = 0.02;

%_______EPI Sessions, for multiple runs with same fieldmap_________________
% Match VDM file to EPI image. This will coregister the field map data to
% the selected EPI for each run/session:
J.matchvdm = 1;

% This will be the name extension followed by an incremented integer for
% run/session specific VDM files:
J.sessname = 'run_';

% Write out distortion corrected EPI image. The image is saved with the 
% prefix 'u'. This is mainly for the quality control of field map correction: 
J.writeunwarped = 1;

% Select an anatomical image for comparison with the distortion corrected EPI:
J.anat = [];

% Match the anatomical image to the distortion corrected EPI:
J.matchanat = 0;

% Adding the EPI sessions/runs images. Either use sbref or first epi_image 
for i=1:numel(epi_files)
    J.session(i).epi ={fullfile(func_dir, epi_files{i})};
end

% Path to the phase image - Change the next two lines if you have a field
% map for each run:
J.phase ={fullfile(fmap_dir, phasediff_img)}; %,'_',num2str(run(1))

% Path to the magnitude image:
J.magnitude =  {fullfile(fmap_dir, magnitude_img)}; %,'_',num2str(run(1))

% Creating the Batch for the SPM:
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj= J;

% Passing the created job to SPM:
spm_jobman('run',matlabbatch);


