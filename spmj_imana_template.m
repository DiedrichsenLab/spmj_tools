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
<<<<<<< Updated upstream
    case 'GLM:make_glm_1'
        
=======
    case 'GLM:make_glm_1'   
        % make the design matrix for the glm
        % models each condition as a separate regressors
        % For conditions with multiple repetitions, one regressor
        % represents all the instances
        % nishimoto_imana('GLM:design1', 'sn', [6])
        
        sn = subj_id;
        hrf_cutoff = Inf;
        prefix = 'r'; % prefix of the preprocessed epi we want to use
        glm = 1;
        vararginoptions(varargin, {'sn', 'hrf_cutoff', 'ses'});
        

        % get the info file that specifies the the tasks and order?
        Dd = dload(fullfile(base_dir, 'task_description.tsv'));
        
        for s = sn
                func_subj_dir = fullfile(base_dir, func_dir,subj_str{s});
 
                % loop through runs within the current sessions
                itaskUni = 0;
                for ses = [1]
                 % create a directory to save the design
                  subj_est_dir = fullfile(base_dir, glm_first_dir,subj_str{s}, sprintf('ses-%02d',ses));
                  dircheck(subj_est_dir)
                  
                  T = []; % task/condition + session + run info
                  J = []; % structure with SPM fields to make the design
                 
                  J.dir            = {subj_est_dir};
                  J.timing.units   = 'secs';
                  J.timing.RT      = 1.3;
                  J.timing.fmri_t  = 16;
                  J.timing.fmri_t0 = 8;
                  
                    % get the list of runs for the current session
                    runs = run_list{ses};
                    for run = 1:2 %length(runs)
                       %V = spm_vol(fullfile(base_dir,func_dir, subj_str{s},sprintf('ses-%02d', ses),sprintf('r%s_run-%02d.nii', subj_str{s}, run)));
                       %numTRs = length(V);
             
                       % get the path to the tsv file
                       tsv_path = fullfile(base_dir, func_dir,subj_str{s});
                       % get the tsvfile for the current run
                       D = dload(fullfile(tsv_path,sprintf('ses-%02d',ses), sprintf('run%d.tsv', run)));
                       
                       % Get the onset and duration of the last sentence
                       lastSentenceOnset = D.onset(end);
                       lastSentenceDuration = D.duration(end);
                       
                       % Convert the end time of the last sentence to TRs
                       endTimeInTRs = ceil((lastSentenceOnset + lastSentenceDuration) / J.timing.RT);


                       % Define scans up to the last sentence's end time
                       N = cell(endTimeInTRs - numDummys, 1);
                       
                       for i = 1:(endTimeInTRs - numDummys)
                           N{i} = fullfile(func_subj_dir, sprintf('ses-%02d', ses), sprintf('%s%s_run-%02d.nii, %d', prefix, subj_str{s}, run, i+numDummys)); % to exclude dummy volumes
                       end % i (image numbers)
                       J.sess(run).scans = N; % scans in the current runs
                        
                       % loop over trials within the current run and build up
                       % the design matrix
                       for ic = 1:length(Dd.task_name)
                           itaskUni = itaskUni+1;
                           % get the indices corresponding to the current
                           % condition.
                           % this line is necessary as there are some
                           % conditions with more than 1 repetition
                           idx = strcmp(D.trial_type, Dd.task_name{ic});
                           fprintf('* %d instances found for condition %s in run %02d\n', sum(idx), Dd.task_name{ic}, run)
                            
                           %
                           % filling in "reginfo"
                           TT.sn        = s;
                           TT.sess      = ses;
                           TT.run       = run;
                           TT.task_name = Dd.task_name(ic);
                           TT.task      = ic;
                           TT.taskUni   = itaskUni;
                           TT.n_rep     = sum(idx);
                            
                           % filling in fields of J (SPM Job)
                           J.sess(run).cond(ic).name = Dd.task_name{ic};
                           J.sess(run).cond(ic).tmod = 0;
                           J.sess(run).cond(ic).orth = 0;
                           J.sess(run).cond(ic).pmod = struct('name', {}, 'param', {}, 'poly', {});
                            
                           % get onset and duration (should be in seconds)
                           onset    = D.onset(idx) - (J.timing.RT*numDummys);
                           fprintf("The onset is %f\n", onset)
                           if onset < 0
                               warning("negative onset found")
                           end
                           duration = D.duration(idx);
                           fprintf("The duration is %f\n", duration);
                            
                           J.sess(run).cond(ic).onset    = onset;
                           J.sess(run).cond(ic).duration = duration;
                            
                           % add the condition info to the reginfo structure
                           T = addstruct(T, TT);
                            
                            
                        end % ic (conditions)
                        
                       J.sess(run).multi     = {''};
                       J.sess(run).regress   = struct('name', {}, 'val', {});
                       J.sess(run).multi_reg = {''};
                       J.sess(run).hpf       = hrf_cutoff; % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
                   end % run (runs of current session)
                
                
               J.fact             = struct('name', {}, 'levels', {});
               J.bases.hrf.derivs = [0 0];
               J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
               J.volt             = 1;
               J.global           = 'None';
               J.mask             = {fullfile(func_subj_dir,'ses-01','rmask_noskull.nii')};
               J.mthresh          = 0.05;
               J.cvi_mask         = {fullfile(func_subj_dir, 'ses-01', 'rmask_gray.nii')};
               J.cvi              =  'fast';
                
               spm_rwls_run_fmri_spec(J);
                
                
               dsave(fullfile(J.dir{1},sprintf('%s_reginfo.tsv', subj_str{s})), T);
               fprintf('- estimates for glm_%d session %d has been saved for %s \n', glm, ses, subj_str{s});
             end % ses (session)
            
            
        end % sn (subject)  
    
    case 'GLM:estimate'      % estimate beta values
        % Example usage: nishimoto_imana('GLM:estimate', 'glm', 1, 'ses', 1, 'sn', 6)
        
        sn       = subj_id; % subject list
        sessions   = [1];       % session number
        
        vararginoptions(varargin, {'sn', 'sessions'})
        
        for s = sn
         
            for ses = sessions
                fprintf('- Doing glm estimation for session %02d %s\n', ses, subj_str{s});
                subj_est_dir = fullfile(base_dir, glm_first_dir,subj_str{s}, sprintf('ses-%02d', ses));         
            
                load(fullfile(subj_est_dir,'SPM.mat'));
                SPM.swd = subj_est_dir;
            
                spm_rwls_spm(SPM);
            end
        end % s (sn),  
         
        
    case 'GLM:T_contrast'    % make T contrasts for each condition
        %%% Calculating contrast images.
        % Example usage: nishimoto_imana('GLM:T_contrast', 'sn', 2, 'glm', 1, 'ses', 1, 'baseline', 'rest')
        
        sn             = subj_id;    % subjects list
        ses            = 1;              % task number
        glm            = 1;              % glm number
        baseline       = 'rest';         % contrast will be calculated against base (available options: 'rest')
        
        vararginoptions(varargin, {'sn', 'glm', 'ses', 'baseline'})
        
        for s = sn
            
            % get the subject id folder name
            fprintf('Contrasts for session %02d %s\n', ses, subj_str{s})
            glm_dir = fullfile(base_dir, glm_first_dir, subj_str{s}, ses_str{ses}); 
            
            cd(glm_dir);
            
            % load the SPM.mat file
            load(fullfile(glm_dir, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = dload(fullfile(glm_dir, sprintf('%s_reginfo.tsv', subj_str{s})));
            
            % t contrast for each condition type
            utask = unique(T.task)';
            idx = 1;
            for ic = utask
                switch baseline
                    case 'myBase' % contrast vs future baseline :)))
                        % put your new contrasts here!
                    case 'rest' % contrast against rest
                        con                          = zeros(1,size(SPM.xX.X,2));
                        con(:,logical((T.task == ic)& (T.n_rep>0))) = 1;
%                         n_rep = length(T.run(T.task == ic));
%                         n_rep_t = T.n_rep(T.task == ic);
%                         name = unique(T.task_name(T.task == ic));
%                         fprintf('- task is %s: \n', name{1});
%                         fprintf('number of reps in all runs = %d\n', n_rep);
%                         fprintf('numberof reps recorded in tsv = %d\n', n_rep_t);
                        con                          = con/abs(sum(con));            
                end % switch base

                % set the name of the contrast
                contrast_name = sprintf('%s-%s', char(unique(T.task_name(T.task == ic))), baseline);
                SPM.xCon(idx) = spm_FcUtil('Set', contrast_name, 'T', 'c', con', SPM.xX.xKXs);
                
                idx = idx + 1;
            end % ic (conditions)
            
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glm_dir, 'SPM_light.mat'), 'SPM')

            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:length(SPM.xCon)
                for n = 1:numel(conName)
                    oldName = fullfile(glm_dir, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName = fullfile(glm_dir, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName, newName);
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn
    
>>>>>>> Stashed changes
end

















