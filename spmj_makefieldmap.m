function spmj_makefieldmap(dataDir, subj_name, run, varargin)
% function spmj_makefieldmap(dataDir, subj_name, run, startTR, varargin)
%   dataDir: data directory of the project (see standard directory
%           structure)
%   run: string for run identifier: 
%           i.e. {'01'  '02','03','04','05','06','07','08'} 
%   subj_name: For directory and filenames, e.g.  's05'
% VARARGINOPTIONS 
%   prefix:  default 'a': Naming is <prefix><subjname>_run<runnumber>.nii
%   image: Number of the image in run to which fieldmap should be aligned
%               (default = 1) 
%   scanType: sub folder in your subject directory
%   subfolderFieldmap: subfolder in the fieldmap directory
%   subfolderRawdata: subfolder in the imaging_data_raw directory
%   rawdataDir:     forces rawdata Directory to different value from the
%               standard naming 

% Tobias Wiestler 2010
% 2010 documentation Joern Diedrichsen
% 06.02.2012 subfolder option replaced with two options 'subfolderFieldmap'
% 'subfolderRawdata'

% 26/April/2012 - Modified by Naveed Ejaz
% Added support for 3D files while keeping backward compatibility for 4D
% files

prefix='a'; 
image=1; 
subfolderRawdata='';
subfolderFieldmap='';
rawdataDir=''; 

use3D=false;

vararginoptions(varargin,{'prefix', 'image', 'subfolderRawdata', 'subfolderFieldmap', 'use3D','rawdataDir'}); 
spm_dir= fileparts(which('spm'));
spmVer=spm('Ver');


% displaying whether using 3D files or not
% disp(['Using 3D: ' num2str(use3D)])

%_______DEFAULTS_________________________________
J.defaults.defaultsval.et = [10 12.43];                                                                 
J.defaults.defaultsval.maskbrain = 1;                                                                   
J.defaults.defaultsval.blipdir = -1;                                                                    
J.defaults.defaultsval.tert = 53.76;                                                                    
J.defaults.defaultsval.epifm = 0;                                                                       
J.defaults.defaultsval.ajm = 0;                                                                         
J.defaults.defaultsval.uflags.method = 'Mark3D';                                                        
J.defaults.defaultsval.uflags.fwhm = 10;                                                                
J.defaults.defaultsval.uflags.pad = 0;                                                                  
J.defaults.defaultsval.uflags.ws = 1;   
switch (spmVer)
    case 'SPM12'
        J.defaults.defaultsval.mflags.template = {fullfile(spm_dir,'canonical','avg152T1.nii')}; 
    case 'SPM8'
        J.defaults.defaultsval.mflags.template = {fullfile(spm_dir,'templates','T1.nii')}; 
end; 
J.defaults.defaultsval.mflags.fwhm = 5;                                                                 
J.defaults.defaultsval.mflags.nerode = 2;                                                               
J.defaults.defaultsval.mflags.ndilate = 4;                                                              
J.defaults.defaultsval.mflags.thresh = 0.5;                                                             
J.defaults.defaultsval.mflags.reg = 0.02;                                                               
J.matchvdm = 1;                                                                                                                                                                              
J.writeunwarped = 1;                                                                                    
J.anat = [];                  
J.matchanat = 0; 


if (isempty(rawdataDir))
    rawdataDir=fullfile(dataDir, 'imaging_data_raw', subj_name, subfolderRawdata); 
end; 

%_______for multiple run with same fieldmap_________________________________
for i=1:numel(run) 
    if use3D
        J.session(i).epi ={fullfile(rawdataDir, [prefix subj_name,'_run',run{i},'_',num2str(image),'.nii'])};
    else
        J.session(i).epi ={fullfile(rawdataDir, [prefix,subj_name,'_run_',run{i},'.nii,',num2str(image)])};    
    end;
end
J.sessname = {'run_RENAMErunNumber'};%{['run_', run{i}]};
%_______change code here if we have 1 fieldmap for each run_________________________________
J.phase ={fullfile(dataDir, 'fieldmaps', subj_name, subfolderFieldmap, [subj_name,'_phase.nii,1'])};          %,'_',num2str(run(1))
J.magnitude =  {fullfile(dataDir, 'fieldmaps', subj_name, subfolderFieldmap, [subj_name,'_magnitude.nii,1'])};  %,'_',num2str(run(1))

matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj= J;

spm_jobman('run',matlabbatch);
