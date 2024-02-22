function varargout = spmj_imana_template(what,varargin)
% Template function for preprocessing of the fMRI data.
% Rename this function to <experiment_code>_imana.m 
% Don't forget to add path the required tools!


%% Directory specification
% Define the base directory of data:
baseDir = somePath;

BIDS_dir = fullfile(baseDir, 'BIDS'); % Raw data post AutoBids conversion
anatomicalDir   = fullfile(baseDir, 'anatomicals'); % Preprocessed anatomical data (LPI + center AC + segemnt)
imagingDirRaw   = fullfile(baseDir, 'imaging_data'); % Preprocessed functional images (realign + coreg + makesamealign)
behaviourDir    = fullfile(baseDir, 'behavioural_data'); % Timing data from the scanner
fs_dir   = fullfile(baseDir,'surfaceFreeSurfer'); %
surfwbDir        = fullfile(baseDir, 'surfaceWB');
glm_first_dir  = fullfile(baseDir, 'GLM_firstlevel');

%% subject info

% Read from participants .tsv file =======================
pinfo = dload(fullfile(baseDir,'participants.tsv'));
subj_id = pinfo.participant_id; 
sn = [1:length(subj_id)];

%% MAIN OPERATION =========================================================

switch(what)
    case 'ANAT:reslice_LPI'
        % Reslices the Anatomical image to the LPI space.
        for s=sn
            fprintf('- Reslice LPI for %s\n', subj_id);
            
             %unzip bids
            anatomical_zip_name = sprintf('%s_ses-01_acq-MP2RAGEpostproc_run-01_T1w.nii.gz',subj_id);
            anatomical_zip = fullfile(BIDS_dir,subj_id,'ses-01/anat/',anatomical_zip_name);
            gunzip(anatomical_zip)
            
            % (1) Reslice anatomical image to set it within LPI co-ordinate frames
            anatomical_unzipped = sprintf('%s_ses-01_acq-MP2RAGEpostproc_run-01_T1w.nii',subj_id);
            source  = fullfile(BIDS_dir,subj_id,'ses-01/anat/',anatomical_unzipped);
            
            % Check if dest_dir exists 
            dest_dir = fullfile(anatomical_dir,subj_id);
            dircheck(dest_dir);
            
            dest    = fullfile(dest_dir,'anatomical_LPI.nii');
            spmj_reslice_LPI(source,'name', dest);
            
            % (2) In the resliced image, set translation to zero
            V               = spm_vol(dest);
            dat             = spm_read_vols(V);
            V.mat(1:3,4)    = [0 0 0];
            spm_write_vol(V,dat);
            display 'Manually retrieve the location of the anterior commissure (x,y,z) before continuing'
        end
    
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
        
    case 'ANAT:segment' % segment the anatomical image
        % check results when done
        SPMhome = fileparts(which('spm.m'));
        J       = []; % spm jobman
        for s = subj_n
            subj_row=getrow(participant_data,participant_data.sn== s );
            subj_id = subj_row.participant_id{1};
            fprintf('- Anatomical segmentation for %s\n', subj_id);
            
            % Get the directory of subjects anatomical
            subj_anat = fullfile(anatomical_dir, subj_id, 'anatomical.nii');
            
          
            
            J.channel.vols     = {subj_anat};
            J.channel.biasreg  = 0.001;
            J.channel.biasfwhm = 60;
            J.channel.write    = [1 0];
            J.tissue(1).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,1')};
            J.tissue(1).ngaus  = 1;
            J.tissue(1).native = [1 0];
            J.tissue(1).warped = [0 0];
            J.tissue(2).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,2')};
            J.tissue(2).ngaus  = 1;
            J.tissue(2).native = [1 0];
            J.tissue(2).warped = [0 0];
            J.tissue(3).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,3')};
            J.tissue(3).ngaus  = 2;
            J.tissue(3).native = [1 0];
            J.tissue(3).warped = [0 0];
            J.tissue(4).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,4')};
            J.tissue(4).ngaus  = 3;
            J.tissue(4).native = [1 0];
            J.tissue(4).warped = [0 0];
            J.tissue(5).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,5')};
            J.tissue(5).ngaus  = 4;
            J.tissue(5).native = [1 0];
            J.tissue(5).warped = [0 0];
            J.tissue(6).tpm    = {fullfile(SPMhome,'tpm/TPM.nii,6')};
            J.tissue(6).ngaus  = 2;
            J.tissue(6).native = [0 0];
            J.tissue(6).warped = [0 0];

            J.warp.mrf     = 1;
            J.warp.cleanup = 1;
            J.warp.reg     = [0 0.001 0.5 0.05 0.2];
            J.warp.affreg  = 'mni';
            J.warp.fwhm    = 0;
            J.warp.samp    = 3;
            J.warp.write   = [1 1];
            matlabbatch{1}.spm.spatial.preproc=J;
            spm_jobman('run',matlabbatch);
        end % s (subject) 
        
    case 'FUNC:realign'          % realign functional images
        % SPM realigns all volumes to the mean volume of first run
        
                
        for s = sn
            spm_jobman('initcfg')
            
            data = {};
                % initialize data cell array which will contain file names for runs/TR images
                func_ses_subj_dir = fullfile(imaging_dir ,subj_id);
                                
                for r = runs
                    % Obtain the number of TRs for the current run
                    for j = 1:numTRs - numDummys
                        data{r}{j,1} = fullfile(func_ses_subj_dir,sprintf('%s_run-%02d.nii,%d', subj_id, r,j));
                    end % j (TRs/images)
                end % r (runs)            
            spmj_realign(data);
            fprintf('- runs realigned for %s  ',subj_id);

        end % s (sn)
        
    case 'FUNC:coreg'            % coregistration with the anatomicals using spm
        % (1) Manually seed the functional/anatomical registration
        % - Select anatomical image and mean functional image to overlay
        % - Manually adjust mean functional image and save the results ("r" will be added as a prefix)
        
        sn       = subj_id;   % list of subjects        
        step     = 'auto';  % first 'manual' then 'auto'
        prefix   = 'r';      % to use the bias corrected version, set it to 'rbb'
        % ===================
        % After the manual registration, the mean functional image will be
        % saved with r as the prefix which will then be used in the
        % automatic registration
        vararginoptions(varargin, {'sn', 'step', 'prefix'});
        spm_jobman('initcfg')
        for s = sn
            % Get the directory of subjects anatomical and functional
            subj_anat_dir = fullfile(anatomical_dir,  subj_id);
            subj_func_dir = fullfile(imaging_dir, subj_id);
                        
            
            % (2) Automatically co-register functional and anatomical images
            J.ref = {fullfile(subj_anat_dir,'anatomical.nii')};
            
            J.source = {fullfile(subj_func_dir, 'ses-01',sprintf('mean%s_run-01.nii', subj_id))};

            
            J.other             = {''};
            J.eoptions.cost_fun = 'nmi';
            J.eoptions.sep      = [4 2];
            J.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            J.eoptions.fwhm     = [7 7];
            
            % 
            matlabbatch{1}.spm.spatial.coreg.estimate=J;
            spm_jobman('run',matlabbatch);
            
        end % s (sn) 
    
    case 'FUNC:make_samealign' % align all the functionals
    % Aligns all functional images to rmean functional image

    for s = sn
        % Get the directory of subjects functional
        subj_func_dir = fullfile(imaging_dir,subj_id);


            fprintf('- make_samealign %s ' ,subj_id)

            % Select image for reference 
            %%% note that functional images are aligned with the first
            %%% run from first session hence, the ref is always rmean<subj>_ses-01_run-01
            P{1} = fullfile(subj_func_dir,'ses-01', sprintf('%smean%s_run-01.nii', prefix,subj_id));

            % Select images to be realigned
            Q = {};
            for r = runs
                for i = 1:numTRs
                    Q{end+1}    = fullfile(subj_func_dir,...
                                           sprintf('ses-%02d', ses),sprintf('r%s%s_run-%02d.nii,%d', prefix, subj_str{s}, r, i));
                end
            end % r(runs)

            spmj_makesamealign_nifti(char(P),char(Q));
    end % s (sn)
    
    case 'FUNC:make_maskImage'   % make mask images (noskull and grey_only)
           % Make maskImage in functional space

           for s = sn
               % Get the directory of subjects anatomical and functional
               subj_anat_dir = fullfile(anatomical_dir, subj_id);
               subj_func_dir = fullfile(imaging_dir,subj_id);

               fprintf('- make mask for %s\n', subj_id);

               nam{1}  = fullfile(subj_func_dir,'ses-01', sprintf('%smean%s_run-01.nii', prefix, subj_str{s}));
               nam{2}  = fullfile(subj_anat_dir, 'c1anatomical.nii');
               nam{3}  = fullfile(subj_anat_dir, 'c2anatomical.nii');
               nam{4}  = fullfile(subj_anat_dir, 'c3anatomical.nii');
               spm_imcalc(nam, fullfile(subj_func_dir,'ses-01','rmask_noskull.nii'), 'i1>0 & (i2+i3+i4)>0.1')
               
               nam     = {};
               nam{1}  = fullfile(subj_func_dir, 'ses-01', sprintf('%smean%s_run-01.nii', prefix, subj_str{s}));
               nam{2}  = fullfile(subj_anat_dir, 'c1anatomical.nii');
               spm_imcalc(nam, fullfile(subj_func_dir,'ses-01','rmask_gray.nii'), 'i1>0 & i2>0.1')

           end
    
    case 'GLM:make_glm_1'    % design glm
        % make the design matrix for the glm
        % models each condition as a separate regressors
        % For conditions with multiple repetitions, one regressor
        % represents all the instances
        
        sn = [1:length(pinfo.participant_id)];
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
                        
                        % Regressors of no interest 
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
    
    case 'SURF:reconall' % Freesurfer reconall routine
        % Calls recon-all, which performs, all of the
        % FreeSurfer cortical reconstruction process
        
        sn   = subj_id; % subject list
        
        vararginoptions(varargin, {'sn'});
        
        % Parent dir of anatomical images    
        for s = sn
            fprintf('- recon-all %s\n', subj_str{s});
                        % Get the directory of subjects anatomical;
            freesurfer_reconall(fs_dir, subj_str{s}, ...
                      fullfile(anatomical_dir, subj_str{s}, 'anatomical.nii'));
        end % s (sn)
        
    case 'SURF:fs2wb'          % Resampling subject from freesurfer fsaverage to fs_LR
        
        sn   = subj_id; % subject list
        res  = 32;          % resolution of the atlas. options are: 32, 164
        hemi = [1, 2];      % list of hemispheres
       
        vararginoptions(varargin, {'sn', 'res', 'hemi'});

        for s = sn 
            % get the subject id folder name
            outDir   = fullfile(baseDir, 'surfaceWB', 'data'); dircheck(outDir);
            surf_resliceFS2WB(subj_str{s}, fs_dir, outDir, 'hemisphere', hemi, 'resolution', sprintf('%dk', res))
        end % s (sn)  

end


%%  =======================Project-specific Cases==================================

switch(what)
    
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
            
        case 'SUIT:isolate_segment'  
        % Segment cerebellum into grey and white matter
        
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});
        
        for s = sn
            fprintf('- Isolate and segment the cerebellum for %s\n', ...
                subj_str{s})
            spm_jobman('initcfg')
            
            % Get the file of subjects anatomical
            anat_subj_dir  = fullfile(anatomical_dir, subj_str{s});
            anat_name = 'anatomical.nii'

            % Define suit folder
            suit_dir = fullfile(baseDir, 'suit/anatomicals',subj_str{s});
            % Create suit folder if it does not exist
            if ~exist(suit_dir, 'dir')
                mkdir (suit_dir)
            end
            
            % Copy anatomical_raw file to suit folder
            source = fullfile(anat_subj_dir, anat_name);
            dest   = fullfile(suit_dir, anat_name);           
            copyfile(source, dest);
            
            % go to subject directory for suit and isolate segment
            suit_isolate_seg({dest}, 'keeptempfiles', 1);
        end % s (sn)
        
        case 'SUIT:normalise_dartel' % SUIT normalization using dartel
        % LAUNCH SPM FMRI BEFORE RUNNING!!!!!
        sn = subj_id; %subjNum
        vararginoptions(varargin, 'sn');
        
        for s = sn
            suit_subj_dir = fullfile(baseDir, 'suit/anatomicals', subj_str{s});
            job.subjND.gray       = {fullfile(suit_subj_dir,'c_anatomical_seg1.nii')};
            job.subjND.white      = {fullfile(suit_subj_dir,'c_anatomical_seg2.nii')};
            job.subjND.isolation  = {fullfile(suit_subj_dir,'c_anatomical_pcereb.nii')};
            suit_normalize_dartel(job);

        end % s (subjects)

        case 'SUIT:save_dartel_def'    
        sn = subj_id; %subjNum
        % Saves the dartel flow field as a deformation file. 
        for s = sn
            cd(fullfile(baseDir,'suit/anatomicals', subj_str{s}));
            anat_name = 'anatomical';
            suit_save_darteldef(anat_name);
        end

        case 'SUIT:reslice'            % Reslice stuff into suit space 
        % run the case with 'anatomical' to check the suit normalization
        % make sure that you reslice into 2mm^3 resolution
        
        sn   = subj_id;
        type = 'con';  % 'betas' or 'con' or 'ResMS' or 'cerebellarGrey' or 'anatomical'
        mask = 'c_anatomical_pcereb'; % 'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask' or 'pcereb'
        glm  = 1;             % glm number. Used for reslicing betas and contrasts 
        
        vararginoptions(varargin, {'sn', 'type', 'mask', 'glm'})
        
        for s = sn
            suit_dir = fullfile(baseDir, 'suit/anatomical',subj_str{s});
            switch type
                case 'anatomical'
                    subj_dir = suit_dir;
                    % Get the name of the anatpmical image
                    files2map = sprintf('%s_T1w_lpi.nii', subj_str{s});
                    
                    job.subj.resample = {sprintf('%s,1', files2map)};
                 case 'betas'
                    glmSubjDir = fullfile(glm_first_dir,sprintf('glm_%d',glm),subj_str{s});
                    images='resbeta_0';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'con'
                    glmSubjDir = fullfile(glm_first_dir,sprintf('glm_%d',glm),subj_str{s});
                    images='con_';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'spmT'
                    glmSubjDir = fullfile(glm_first_dir,sprintf('glm_%d',glm),subj_str{s});
                    images='spmT_';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
            end
            job.subj.affineTr = {fullfile(baseDir,'suit','anatomicals',subj_str{s},'Affine_c_anatomical_seg1.mat')};
            job.subj.flowfield= {fullfile(baseDir,'suit','anatomicals',subj_str{s},'u_a_c_anatomical_seg1.nii')};
            job.subj.resample = {source.name};
            job.subj.mask     = {fullfile(baseDir,'suit','anatomicals',subj_str{s},sprintf('%s.nii',mask))};
            job.vox           = [1 1 1];
            % Replace Nans with zeros to avoid big holes in the the data 
            for i=1:length(source)
                V=spm_vol(source(i).name); 
                X=spm_read_vols(V); 
                X(isnan(X))=0; 
                spm_write_vol(V,X); 
            end; 
            suit_reslice_dartel(job);
            
            source=fullfile(glmSubjDir,'*wd*');
            destination=fullfile(baseDir,'suit',sprintf('glm_%d',glm),subj_str{s});
            movefile(source,destination);

            fprintf('%s have been resliced into suit space \n',type)
        end







