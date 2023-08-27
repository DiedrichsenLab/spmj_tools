function [A,hdr,mf,data_type] = imgload (filename,varargin)
% [A,hdr]=imgload (name, varargin) 
% filename: either with or without explicit *.img ending 
%           looks for the corresponding *.hdr file
% uses spm_read_hdr to asess the header information
%   this information is then used to put the image into a data matrix 
%   'parametername',new_value 
%   can be used in varargin to change any of the settings 
if (strcmp(filename(end-3:end),'.img'))
    filename=filename(1:end-4);
end;
[hdr,otherendian] = spm_read_hdr([filename '.hdr']);
if (isempty(hdr))
    error(sprintf('file not found: %s',[filename '.hdr'])); 
end;
if (otherendian)
    mf='ieee-be';
else
    mf='ieee-le';
end;    
%data_type=deblank(hdr.hk.data_type);
%if isempty(data_type)
switch (hdr.dime.datatype)
    case 0
        data_type='int16'; % ???? doubtful
    case 2
        data_type='int8';
    case 4
        data_type='int16';
    case 16
        data_type='float32';        
    case 64
        data_type='float64';
end;
    %end;
fid=fopen([filename '.img'],'r',mf);
[brain,count]=fread(fid,inf,data_type);
fclose(fid);
if (hdr.dime.dim(1)==3)
    A=reshape(brain,hdr.dime.dim(2),hdr.dime.dim(3),hdr.dime.dim(4));
elseif (hdr.dime.dim(1)==4)
    if (hdr.dime.dim(5)==1)
        A=reshape(brain,hdr.dime.dim(2),hdr.dime.dim(3),hdr.dime.dim(4));
    else
        fprintf('Four Dimensions\n');
    end;    
end;    
