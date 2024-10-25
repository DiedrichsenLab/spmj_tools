function change_nifti_numformat(infile,outfile, varargin); 
% Changes the number format of a nifti-file. 
%    Args:
%        infile (str): file name of input file 
%        outfile (str): file name output file 
%    Varargin:
%        new_numformat (str): New number format. Defaults to 'uint16'.
%        typecast_data (bool): If true, typecasts the data as well. 
%     Note:
%        typecast_data changes the data format and can lead to and wrap-around of values, as the data is typecasted as well. 
%        Do only when correcting a previous mistake in processing. 
new_numformat = 'uint16';  
typecast_data = 1; 
V = spm_vol(infile);
if (typecast_data)
    V.dt = [spm_type(new_numformat) 0]; 
end; 
X = spm_read_vols(V,0);
V.fname = outfile; 
V.dt = [spm_type(new_numformat) 0]; 
spm_write_vol(V,X); 