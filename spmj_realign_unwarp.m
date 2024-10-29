function spmj_realign_unwarp(subj_name, run_names, varargin)
% spmj_realign_unwarp(subj_name, run_names)
% INPUT: 
%   subj_name:  Name for subdirectory and file-prefix (i.e. 's05')
%   run_names:  Cell array of identifiers for the run 
%               {'run-01' 'run-02','run-03','run-04'} 
% VARARGINOPTIONS: 
%   'sess_names':       optional session name if multiple scannning sesss=subfoldes are used 
%   'scanType'          sub folder in your subject directory
%   'fmap_dir'          folder in the fieldmap directory
%   'rawdata_dir'       folder in the unaligned imaging_data_raw directory
%   'base_dir':         Root directory for the imaging structure (needs directories 
%                       imaging_data_raw and fieldmaps 
%   'raw_name':         Different name for file ('run', vs. 'sbref')
%                       defaults to 'run' 
% EXAMPLE USAGE: 
% A. Single session experiment:
% leave sesssion names empty, run_names is a single cell array of images
% B.Multi session experiment:
% sess_names = {'ses-01','ses-02}
% run_names = {{'run-01' 'run-02','run-03','run-04' },{'run-01' 'run-02'}}
raw_name = 'run'; 
rawdata_dir = ''; 
fmap_dir = ''; 
base_dir = ''; 
sess_names = {}; 
rtm = 0;    % register to mean or first volume. default register to first volume of the first run.

vararginoptions(varargin,{'prefix', 'fmap_dir','rawdata_dir','base_dir','target_dir','rtm','raw_name','sess_names'}); 

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

if isempty(rawdata_dir)
    rawdata_dir = fullfile(base_dir,'imaging_data_raw');
end
if isempty(fmap_dir)
    fmap_dir = fullfile(base_dir,'fieldmaps');
end 

%_______images and fieldmap definition_________________________________
% Flat hierarchy one imaging session 
if (isempty(sess_names))
    for j=1:numel(run_names)
        file_name = fullfile(rawdata_dir,subj_name,sprintf('%s_%s.nii',subj_name,run_names{s}{j}));
        if (V.dat.dim==3)
            scans{1}= file_name;
        else 
            imageNumber=1:V.dat.dim(4); 
            for i= 1:numel(imageNumber)
                scans{i}= sprintf('%s,%d',file_name,i);
            end
        end
        J.data(j).scans = scans';
        clear scans; 
        J.data(j).pmscan = {fullfile(dataDir, 'fieldmaps',subj_name,subfolderFieldmap,['vdm5_sc',subj_name,'_phase_',run_names{j},'.nii,1'])};
    end 
else  % Nested hierarchy with different subfolders for different sessions 
    indx = 1; 
    for s=1:numel(sess_names)
        for j=1:numel(run_names{s})
            file_name = fullfile(rawdata_dir,subj_name,sess_names{s},sprintf('%s_%s.nii',subj_name,run_names{s}{j}));
            V = nifti(file_name);
            if (length(V.dat.dim)==3) % 3dimage (likely sbref)
                scans{1}= file_name;
            else 
                imageNumber=1:V.dat.dim(4); 
                for i= 1:numel(imageNumber)
                    scans{i}= sprintf('%s,%d',file_name,i);
                end
            end
            J.data(indx).scans = scans';
            J.data(indx).pmscan = {fullfile(fmap_dir,subj_name,sess_names{s},sprintf('vdm5_sc%s_%s_phasediff_run_%d.nii',subj_name,sess_names{s},j))};
            indx = indx+1; 
            clear scans; 
        end
    end
end

matlabbatch{1}.spm.spatial.realignunwarp= J;
spm_jobman('run',matlabbatch);



