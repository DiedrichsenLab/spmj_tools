function spmj_write_hdr(fname,hdr,swapped);
% function spmj_write_hdr(fname,hdr,swapped);
% Does the opposite of spm_read_hdr
% Just writes header information to disk, without 
% chaning any of the fields as SPM does. 
% However: all fields have to be given 
if spm_platform('bigend')
    if (swapped)
        mach='ieee-le';
    else 
        mach='ieee-be';
    end;
else
    if (swapped)
        mach='ieee-be';
    else 
        mach='ieee-le';
    end;
end;
fid           = fopen(fname,'w',mach);
if (fid == -1),
	error(['Error opening ' fname '. Check that you have write permission.']);
end;

write_hk(fid,hdr.hk);
write_dime(fid,hdr.dime);
write_hist(fid,hdr.hist);
fclose(fid);
%_______________________________________________________________________
%_______________________________________________________________________
function write_hk(fid,hk)
% write (struct) header_key
%-----------------------------------------------------------------------
fseek(fid,0,'bof');
fwrite(fid,hk.sizeof_hdr,	'int32');
fwrite(fid,hk.data_type,	'char' );
fwrite(fid,hk.db_name,		'char' );
fwrite(fid,hk.extents,		'int32');
fwrite(fid,hk.session_error,'int16');
fwrite(fid,hk.regular,		'char' );
if fwrite(fid,hk.hkey_un0,	'char' )~= 1,
	error(['Error writing '  fopen(fid) '. Check your disk space.']);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function write_dime(fid,dime)
% write (struct) image_dimension
%-----------------------------------------------------------------------
fseek(fid,40,'bof');
fwrite(fid,dime.dim,		'int16');
fwrite(fid,dime.vox_units,	'uchar' );
fwrite(fid,dime.cal_units,	'uchar' );
fwrite(fid,dime.unused1,	'int16' );
fwrite(fid,dime.datatype,	'int16');
fwrite(fid,dime.bitpix,		'int16');
fwrite(fid,dime.dim_un0,	'int16');
fwrite(fid,dime.pixdim,		'float');
fwrite(fid,dime.vox_offset,	'float');
fwrite(fid,dime.funused1,	'float');
fwrite(fid,dime.funused2,	'float');
fwrite(fid,dime.funused2,	'float');
fwrite(fid,dime.cal_max,	'float');
fwrite(fid,dime.cal_min,	'float');
fwrite(fid,dime.compressed,	'int32');
fwrite(fid,dime.verified,	'int32');
fwrite(fid,dime.glmax,		'int32');
if fwrite(fid,dime.glmin,		'int32')~=1,
	error(['Error writing '  fopen(fid) '. Check your disk space.']);
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function write_hist(fid,hist)
% write (struct) data_history
%-----------------------------------------------------------------------
fseek(fid,148,'bof');
fwrite(fid,hist.descrip,	'uchar');
fwrite(fid,hist.aux_file,	'uchar');
fwrite(fid,hist.orient,		'uchar');
fwrite(fid,hist.origin,		'int16');
fwrite(fid,hist.generated,	'uchar');
fwrite(fid,hist.scannum,	'uchar');
fwrite(fid,hist.patient_id,	'uchar');
fwrite(fid,hist.exp_date,	'uchar');
fwrite(fid,hist.exp_time,	'uchar');
fwrite(fid,hist.hist_un0,	'uchar');
fwrite(fid,hist.views,		'int32');
fwrite(fid,hist.vols_added,	'int32');
fwrite(fid,hist.start_field,'int32');
fwrite(fid,hist.field_skip,	'int32');
fwrite(fid,hist.omax,		'int32');
fwrite(fid,hist.omin,		'int32');
fwrite(fid,hist.smax,		'int32');
if fwrite(fid,hist.smin,		'int32')~=1,
	error(['Error writing '  fopen(fid) '. Check your disk space.']);
end;
return;