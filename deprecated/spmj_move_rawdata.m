function SPM=spmj_move_rawdata(SPM,rawdata_dir)
% function SPM=spmj_move_rawdata(SPM,rawdata_dir)
% Updates the field SPM.xY.P and SPM.xY.VY to access the raw data under a
% different directory 
% INPUT:  
%   SPM:              SPM or filename. If Filename is given it's loaded and saved
%   rawdata_dir:      New directory for the raw data 
% OUTPUT: 
%   SPM:              Updated SPM structure (is saved automatically, if
%                     filename
% j.diedrichsen@ucl.ac.uk 
if (ischar(SPM))
    SPM_name=SPM; 
    load(SPM_name); 
end; 

% Update file names 
Rawcell={};
for i=1:size(SPM.xY.P,1)
    [d,name,ext,num]=spm_fileparts(SPM.xY.P(i,:));
    if (isempty(d)) % No directory - we think this is maybe from a windows platform (\)
        indx=find(name=='\'); 
        name=name(indx(end)+1:end); 
    end; 
    Rawcell{i}=fullfile(rawdata_dir,[name ext num]);
end;

% Place and update volumes 
SPM.xY.P  = char(Rawcell);
SPM.xY.VY = spm_vol(SPM.xY.P); 

% Save SPM if necessary 
if (exist('SPM_name','var'))
    save SPM_name SPM; 
end; 
