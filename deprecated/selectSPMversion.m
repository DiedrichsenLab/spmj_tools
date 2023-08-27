function val = selectSPMversion(version);
%% function val = selectSPMversion(version);
% Select spm version and change path accordingly.
% 
% Usage: selectSPMversion('spm12');
% 
%        if (selectSPMversion('spm12'))
%           % do something
%        end
% 
% a-yokoi (at.yokoi.work@gmail.com)
val = true;

% 1. Get current version of spm in the search path
spmDir      = which('spm.m');
tmp         = strsplit(spmDir,filesep);
currVersion = tmp{end-1}; % =spm('Ver')

% 2. Check if specified spm version exists in the computer 
%   (assuming they are all in the same directory)
cd(fullfile(filesep,tmp{1:end-2}));
if ~exist(version,'dir')
    warning('%s doesn''t exist!!!',version);
    val = false;
    return;
end
newVersion = version;

% 3. Remove older version from the search path
currPath    = strsplit(path,':'); % matlab seems to use : to separate path
spmPath     = strfind(currPath,currVersion);
spmPath     = currPath(~cellfun(@isempty,spmPath));
rmpath(spmPath{:});

% 4. Add new version to the search path (doesn't save path)
cd(spmDir(1:end-5));
cd ../
baseDir = cd;
addpath(genpath(fullfile(baseDir,newVersion))); % don't need to add all sub-folders?
% addpath(fullfile(baseDir,newVersion));

end