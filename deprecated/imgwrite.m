% function imgwrite (name,A,format)
% writes an image-file to disk
function imgwrite (name,A,format)
  if(nargin<3)
      format='int16';
  end;
  fid=fopen(name,'w');
  fwrite(fid,A,format);
  fclose(fid);

