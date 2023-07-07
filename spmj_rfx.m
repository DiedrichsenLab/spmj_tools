function SPM=spmj_rfx(filename,subjnames,varargin) 


% Digest input parameters 
PARAM=struct(varargin)
if exist('PARAM', 'var') & ~isempty(PARAM)
    fnames = fieldnames(PARAM);
    for fnum = 1:length(fnames);
        eval([fnames{fnum} ' = getfield(PARAM, fnames{fnum});']);
    end;
end;

