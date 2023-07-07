function spmj_regress_parameter(dataDir, images, covariate, covariate_names) 
%TW2011
J.dir = {dataDir};                             

J.des.mreg.scans = images;                                                                            

for i= 1:size(covariate, 2)
    J.des.mreg.mcov(i).c = covariate(:,i);
    J.des.mreg.mcov(i).cname = covariate_names{i};
    J.des.mreg.mcov(i).iCC = 5;
end

J.des.mreg.incint = 0;                                                                            
J.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});                                      
J.masking.tm.tm_none = 1;                                                                         
J.masking.im = 1;                                                                                 
J.masking.em = {''};                                                                              
J.globalc.g_omit = 1;                                                                             
J.globalm.gmsca.gmsca_no = 1;                                                                     
J.globalm.glonorm = 1; 

matlabbatch{1}.spm.stats.factorial_design= J;

spm_jobman('run',matlabbatch);
