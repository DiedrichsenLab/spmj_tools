function spmj_realign_estimate(dataDir, subj_name, run, startTR, endTR, varargin)

%toDo docu
% Tobias Wiestler 2010

prefix= 'a';
scanType=[];
vararginoptions(varargin,{'prefix', 'scanType'}); 
 
imageNumber= startTR:endTR;

for j=1:numel(run)
    for i= 1:numel(imageNumber)
        session_run{i}= fullfile(dataDir, 'imaging_data_raw',subj_name,scanType, [prefix,subj_name,'_run',run{j},'.nii,',num2str(imageNumber(i))]);                              
    end
    data{j}= session_run';
end 
J.data = data'; 

J.eoptions.quality = 0.9;                                                                               
J.eoptions.sep = 4;                                                                                     
J.eoptions.fwhm = 5;                                                                                    
J.eoptions.rtm = 1;                                                                                     
J.eoptions.interp = 2;                                                                                  
J.eoptions.wrap = [0 0 0];                                                                              
J.eoptions.weight = {''};                                                                               
J.roptions.which = [2 1];                                                                               
J.roptions.interp = 4;                                                                                  
J.roptions.wrap = [0 0 0];                                                                              
J.roptions.mask = 1;                                                                                    
J.roptions.prefix = 'r'; 


matlabbatch{1}.spm.spatial.realign.estwrite=J;
spm_jobman('run',matlabbatch);
