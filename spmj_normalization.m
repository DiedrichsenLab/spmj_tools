function spmj_normalization(greymatter_image)
%TO CHECK
J.subj.source = {[greymatter_image,',1']};
J.subj.wtsrc = '';         
spmLocation= which('SPM');
J.eoptions.template = {fullfile(spmLocation(1:end-6),'templates','T1.nii,1')};                       
J.eoptions.weight = '';                                                                       
J.eoptions.smosrc = 8;                                                                        
J.eoptions.smoref = 0;                                                                        
J.eoptions.regtype = 'mni';                                                                   
J.eoptions.cutoff = 25;                                                                       
J.eoptions.nits = 16;                                                                         
J.eoptions.reg = 1;  

matlabbatch{1}.spm.spatial.normalise.est=J;
spm_jobman('run',matlabbatch);