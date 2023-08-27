function spmj_write_normalization(deformation_image, images)
%TO CHECK
J.subj.matname = {deformation_image};
J.subj.resample = cellstr(images);   
J.roptions.preserve = 0;                                                                          
J.roptions.bb = [-78 -112 -50                                                                     
                  78   76  85];                                                                       
J.roptions.vox = [2 2 2];                                                                         
J.roptions.interp = 1;                                                                            
J.roptions.wrap = [0 0 0];                                                                        
J.roptions.prefix = 'w';    

matlabbatch{1}.spm.spatial.normalise.write=J;
spm_jobman('run',matlabbatch);