function spmj_realign_unwarp_sess(dataDir, subj_name, run, nTR, varargin)
% spmj_realign_unwarp(dataDir, subj_name, run, startTR, endTR)
% INPUT: 
%   dataDir:    Root directory for the imaging structure (needs directories 
%               imaging_data_raw and fieldmaos 
%   subj_name:  Name for subdirectory and file-prefix (i.e. 's05') 
%   run:        Cell array of identifiers for the run: 
%                   e.g.: {'01' '02','03','04','05','06','07','08'} 
%               if multiple sessions pass cell array with numSession cells: 
%                   e.g.: {{'01' '02','03','04','05','06','07','08'},...
%                          {'09' '10','11'}};
% VARARGINOPTIONS: 
%   'prefix'            prefix for run name (default 'a'); 
%   'scanType'          sub folder in your subject directory
%   'subfolderFieldmap' subfolder in the fieldmap directory: for multiple
%                       sessions, pass NumSession directories  
%   'subfolderRawdata'  subfolder in the imaging_data_raw directory: for
%                       multiple sessions, pass NumSession directory
%   'rawdataDir'        overwrites standard naming and forces routine to
%                       use this folder for location of raw data 
% Tobias Wiestler & Joern Diedrichsen & Naveed Ejaz 
% 06/02/2012 subfolder option replaced with two options 'subfolderFieldmap' 'subfolderRawdata'
% 23/10/2012 added rawdataDir to be able to overwrite the standard naming convention
% 
prefix= 'a';
subfolderRawdata='';
subfolderFieldmap='';
use3D=false;
rawdataDir=''; 

vararginoptions(varargin,{'prefix', 'subfolderRawdata', 'subfolderFieldmap','use3D','rawdataDir'}); 

% displaying whether using 3D files or not
disp(['Using 3D: ' num2str(use3D)]);


%_______DEFAULTS_________________________________
J.eoptions.quality = 0.9;
J.eoptions.sep = 2;%4;                                                                                   
J.eoptions.fwhm = 5;                                                                                  
J.eoptions.rtm = 0; %why zero and not one                                                                                  
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
J.uwroptions.uwwhich = [2 1]; %[2 1] with mean image without [2 0]
J.uwroptions.rinterp = 4;                                                                             
J.uwroptions.wrap = [0 1 0];  %  wrap-around in the [x y z] direction during the reslicing (wrap)                                                                                
J.uwroptions.mask = 1;                                                                                
J.uwroptions.prefix = 'u'; 

imageNumber= 1:nTR;

if (isempty(rawdataDir))
%     rawdataDir=fullfile(dataDir, 'imaging_data_raw',subj_name,subfolderRawdata); 
    rawdataDir=fullfile(dataDir, 'imaging_data_raw',subj_name); 
end; 

%_______images and fieldmap definition_________________________________
sCount=1;
for s = 1:numel(run)
    for j=1:numel(run{s})
        for i= 1:numel(imageNumber)
            if use3D
                scans{i}= fullfile(rawdataDir,subfolderRawdata{s}, [prefix subj_name,'_run',run{s}{j},'_',num2str(imageNumber(i)),'.nii']);
            else
                scans{i}= fullfile(rawdataDir,subfolderRawdata{s}, [prefix,subj_name,'_run',run{s}{j},'.nii,',num2str(imageNumber(i))]);
            end;

            %['/media/SECOND/SequenceLearning/test_fieldmap/dicom/a', subj_name,'_run',run{j},'_',num2str(imageNumber(i)),'.nii,1'];
            %                                      
        end
        J.data(sCount).scans = scans;
        J.data(sCount).pmscan = {fullfile(dataDir, 'fieldmaps',subj_name,subfolderFieldmap{s},['vdm5_sc',subj_name,'_phase_session',num2str(j),'.nii,1'])};
        sCount=sCount+1;
        %{['/media/SECOND/SequenceLearning/test_fieldmap/dicom/',['vdm5_sc',subj_name,'_phase_1_session',num2str(j),'.nii,1']]};
    end 
end;

matlabbatch{1}.spm.spatial.realignunwarp= J;
spm_jobman('run',matlabbatch);