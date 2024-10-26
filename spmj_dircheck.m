function spmj_dircheck(dir)
% Checking for directories and creating if not exist
if ~exist(dir,'dir')
    warning('%s doesn''t exist. Creating one now.\n',dir);
    mkdir(dir);
end

