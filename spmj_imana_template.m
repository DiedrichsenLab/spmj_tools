function varargout = spmj_imana_template(what,varargin)
% Template function for preprocessing of the fMRI data.
% Rename this function to <experiment_code>_imana.m 
% Don't forget to add path the required tools!


%% Directory specification
% Define the base directory of data:
baseDir = fullfile(somePath,'data/');

fieldmapsDir    = fullfile(baseDir, 'fieldmaps');
behaviourDir    = fullfile(baseDir, 'behavioural_data');
analyzeDir 		= fullfile(baseDir, 'analyze');
anatomicalDir   = fullfile(baseDir, 'anatomicals');
imagingDirRaw   = fullfile(baseDir, 'imaging_data_raw');
imagingDir      = fullfile(baseDir, 'imaging_data');
glmDir          = fullfile(baseDir, 'glm_firstlevel_1');
freesurferDir   = fullfile(baseDir, 'surfaceFreesurfer');
caretDir        = fullfile(baseDir, 'surfaceCaret');
regDir          = fullfile(baseDir, 'RegionOfInterest');
physioDir       = fullfile(baseDir,'physio_data');

%% subject info

% Read from participants .tsv file =======================
pinfo = datload('data/participants.tsv');


%% MAIN OPERATION =========================================================

switch(what)
    case 'ANAT:reslice_LPI'
        % Reslices the Anatomical image to the LPI space. This is not
        % needed for the CFMM scans since they already come in LPI space.
    
    case 'ANAT:recenter_AC'
        % Description:
        % Recenters the anatomical data to the Anterior Commissure
        % coordiantes. Doing that, the [0,0,0] coordiante of subject's
        % anatomical image will be the Anterior Commissure.

        % You should manually find the voxel coordinates of AC 
        % for each from their anatomical scans and add it to the
        % participants.tsv file under the loc_AC column.

        % This function runs for all subjects and sessions.
        

        % location of AC as an Nx3 array, N being number of subjs:
        loc_AC = pinfo.loc_AC;

        % looping through subjects:
        for sn = 1:length(subj_runs)
            % path to the raw anatomical .nii file:
            img_path = fullfile(anatomicalDir,subj_name{sn},strcat(subj_name{sn},'_anatomical_raw.nii'));
            
            % Get header information for the image:
            V = spm_vol(img_path);

            % Reads the image:
            dat = spm_read_vols(V);

            % V.mat is the 4x4 affine transform from index 
            % to real world coordinates. So V.mat(1:3,4) is the 
            % translation vector:
            oldOrig = V.mat(1:3,4);

            % changing the translation vector to put AC at [0,0,0]:
            V.mat(1:3,4) = oldOrig+loc_AC(sn,:)';

            % writing and saving the volume:
            spm_write_vol(V,dat);
            fprintf('recenter AC done for %s \n',subj_name{sn})
        end

        
    case 'FMAP:makefieldmap'
        % Description:
        % Generates VDM files from the presubtracted phase & magnitude
        % images acquired from the field map sequence. Also, just as a
        % quality control this function creates unwarped EPIs from the
        % functional data with the prefix 'u' for each run.

        % This function runs for all subjects and sessions.
        % creating the run for each subject and adding to subj_runs:
        run_file = participants_info.run_sess;
        nRun = 10;
        subj_runs = cell(size(run_file,1),1);
        for i = 1:size(run_file,1)
            [~,ia,~] = unique(run_file(i,:));
            run = {cellfun(@num2str, num2cell([1:ia(2)-1]), 'UniformOutput', false), cellfun(@num2str, num2cell([ia(2):nRun]), 'UniformOutput', false)};
            subj_runs{i} = run;
        end

        % Prefix of the functional files:
        prefixepi  = '';

        % Prefix of the fieldmap files:
        prefixfieldmap  = '';

        % echo times of the gradient eho sequence:
        et1 = 4.92;
        et2 = 7.38;

        % total EPI readout time = = echo spacing (in ms) * base resolution 
        % (also knows as number of echos). If you use GRAPPA acceleration, 
        % you need to divide the total number of echos by two:
        tert = 90 * 0.7 / 2;
        
        % looping through subjects:
        for sn = 1:length(subj_runs)
            % looping through sessions, length(run) = 2 = num sessions:
            for sess = 1:length(run)
                subfolderRawdata = sprintf('sess%d',sess);
                subfolderFieldmap = sprintf('sess%d',sess);
                % function to create the makefieldmap job and passing it to the SPM
                % job manager:
                spmj_makefieldmap(baseDir,subj_name{sn},subj_runs{sn}{sess}, ...
                          'et1', et1, ...
                          'et2', et2, ...
                          'tert', tert, ...
                          'prefix',prefixepi, ...
                          'subfolderRawdata',subfolderRawdata, ...
                          'subfolderFieldmap',subfolderFieldmap);
            end
        end
    case 'GLM:make_glm_1'
        
end

















