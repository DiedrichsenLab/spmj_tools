function spmj_realign_unwarp(dataDir, subj_name, run, startTR, endTR, varargin)
% spmj_realign_unwarp(dataDir, subj_name, run, startTR, endTR)
% INPUT: 
%   dataDir:    Root directory for the imaging structure (needs directories 
%               imaging_data_raw and fieldmaos 
%   subj_name:  Name for subdirectory and file-prefix (i.e. 's05') 
%   run:        Cell array of identifiers for the run 
%               {'01' '02','03','04','05','06','07','08'} 
%   startTR:    First image to align 
%   endTR:      Last image to align (if INF, it will use all available)
%               Array specifying a last image for each run or single value 
%               for all runs 
% VARARGINOPTIONS: 
%   'prefix'            prefix for run name (default 'a'); 
%   'scanType'          sub folder in your subject directory
%   'subfolderFieldmap' subfolder in the fieldmap directory
%   'subfolderRawdata'  subfolder in the imaging_data_raw directory
%   'rawdataDir'        overwrites standard naming and forces routine to
%                       use this folder for location of raw data 
% Tobias Wiestler & Joern Diedrichsen
% 06/02/2012 subfolder option replaced with two options 'subfolderFieldmap' 'subfolderRawdata'
% 23/10/2012 added rawdataDir to be able to overwrite the standard naming convention
% Sungshin Kim
% 24/07/05 endTR option specifying a last image for each run
%
if startTR>1
    error('Set startTR as 1');
end

prefix= 'a';
subfolderRawdata='';
subfolderFieldmap='';
use3D=false;
rawdataDir=''; 
rtm = 0;    % register to mean or first volume. default register to first volume of the first run.

vararginoptions(varargin,{'prefix', 'subfolderRawdata', 'subfolderFieldmap','use3D','rawdataDir','rtm'}); 


%_______DEFAULTS_________________________________
J.eoptions.quality = 0.9;
J.eoptions.sep = 2;%4;                                                                                   
J.eoptions.fwhm = 5;                                                                                  
J.eoptions.rtm = rtm; % register to mean -> 0: register to first volume of each run. 1: register to mean image of each run.                                                                          
J.eoptions.einterp = 2;                                                                               
J.eoptions.ewrap = [0 1 0];     %  wrap-around in the [x y z] direction during the estimation (ewrap)  wrap of the front of the head to the back of the head                                                                       
J.eoptions.weight = {''};                                                                             
J.uweoptions.basfcn = [12 12];                                                                        
J.uweoptions.regorder = 1;                                                                            
J.uweoptions.lambda = 100000;                                                                         
J.uweoptions.jm = 0;                                                                                  
J.uweoptions.fot = [4 5];                                                                             
J.uweoptions.sot = [1];                                                                                 
J.uweoptions.uwfwhm = 4;                                                                              
J.uweoptions.rem = 1;                                                                                 
J.uweoptions.noi = 5;                                                                                 
J.uweoptions.expround = 'Average';                                                                    
J.uwroptions.uwwhich = [2 1]; %[2 1]: with mean image. [2 0]: without 
J.uwroptions.rinterp = 4;                                                                            
J.uwroptions.wrap = [0 1 0];  %  wrap-around in the [x y z] direction during the reslicing (wrap)                                                                                
J.uwroptions.mask = 1;                                                                                
J.uwroptions.prefix = 'u'; 


if (isempty(rawdataDir))
    rawdataDir=fullfile(dataDir, 'imaging_data_raw',subj_name,subfolderRawdata); 
end; 

if numel(endTR)==1 & ~isinf(endTR)
    endTR = endTR*ones(1,numel(run));
end
%_______images and fieldmap definition_________________________________
for j=1:numel(run)
    if (isinf(endTR))  % All avaialble images: only works with 4d-nifits right now 
        V = nifti(fullfile(rawdataDir,[prefix,subj_name,'_run_',run{j},'.nii'])); 
        imageNumber=startTR:V.dat.dim(4); 
    else
        imageNumber= startTR:endTR(j);
    end; 
    for i= 1:numel(imageNumber)
        if use3D
            scans{i}= fullfile(rawdataDir, [prefix subj_name,'_run',run{j},'_',num2str(imageNumber(i)),'.nii']);
        else
            scans{i}= fullfile(rawdataDir, [prefix,subj_name,'_run_',run{j},'.nii,',num2str(imageNumber(i))]);
        end;
    end;
    J.data(j).scans = scans';
    clear scans; %% added by SKim on 20240704
    J.data(j).pmscan = {fullfile(dataDir, 'fieldmaps',subj_name,subfolderFieldmap,['vdm5_sc',subj_name,'_phase_run_',num2str(j),'.nii,1'])};
end

matlabbatch{1}.spm.spatial.realignunwarp= J;
spm_jobman('run',matlabbatch);