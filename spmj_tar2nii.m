function spmj_tar2nii (tar_name, nii_name, varargin)
startTR= 1; 
use3D = 0;
vararginoptions(varargin,{'startTR', 'use3D'}); 

%toDo Docu
% TW 2010

%________________________________________________________
% maybe check if there is one tmp dir already 
fprintf(1,'untar \n'); 
fnames=untar(tar_name);
%________________________________________________________
% get all img 
fprintf(1,'read header and produce nii\n'); 
hdr= spm_dicom_headers(char(fnames));
nii= spm_dicom_convert(hdr,'all','flat','nii');
%________________________________________________________
% sort the img because the cells are not in the correct order!
ordered_nii_names= sort(nii.files);    
%________________________________________________________
if ~use3D
    fprintf(1,'make 4d nii'); 
    spm_file_merge(ordered_nii_names(startTR:end),nii_name);
else 
    fprintf(1,'rename 3d nii');
    j=1;
    for i=startTR:numel(ordered_nii_names)
        [~,name] = fileparts(nii_name);
        movefile(ordered_nii_names{i}, [name, '_', num2str(j), '.nii']);
        j=j+1;
    end
end
    
%________________________________________________________
%delete
fprintf(1,'delete ima\n');
if ~use3D
   todelete=cell2struct([fnames,nii.files'] , 'name');
   delete(todelete.name);
else
   todelete = [fnames'; ordered_nii_names(1:startTR-1)];
   delete(todelete{:});
end


