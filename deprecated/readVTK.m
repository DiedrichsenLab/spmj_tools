function V = readVTK(vtkfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: V = readVTK(vtkfile)
%
%   V:       The matrix to be stored
%   vtkfile: The filename
% 
% Erik Vidholm 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = 0;

% open file
fid = fopen(vtkfile,'r','ieee-be');

if( fid == -1 )
  return
end

fgetl(fid); % # vtk DataFile Version x.x
fgetl(fid); % comments
fgetl(fid); % BINARY
fgetl(fid); % DATASET STRUCTURED_POINTS

s = fgetl(fid); % DIMENSIONS NX NY NZ
sz = sscanf(s, '%*s%d%d%d').'

fgetl(fid); % ORIGIN OX OY OZ
fgetl(fid); % SPACING SX SY SZ
fgetl(fid); % POINT_DATA NXNYNZ

s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
svstr = sscanf(s, '%s', 1)
dtstr = sscanf(s, '%*s%*s%s')

if( strcmp(svstr,'SCALARS') > 0 )
  fgetl(fid); % the lookup table
  if( strcmp(dtstr,'unsigned_char') > 0 ) 
    % read data
    V = fread(fid,prod(sz),'*uint8');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'unsigned_short') > 0 )
    % read data
    V = fread(fid,prod(sz),'*uint16');
    V = reshape(V,sz);
  elseif( strcmp(dtstr,'float') > 0 )
    % read data
    V = fread(fid,prod(sz),'*single');
    V = reshape(V,sz);
  end
  
elseif( strcmp(svstr,'VECTORS') > 0 )
  if( strcmp(dtstr,'float') > 0 ) 
    % read data
    V = fread(fid,3*prod(sz),'*single');
    V = reshape(V,[3 sz]);
    V = permute(V,[2 3 4 1]);
    elseif( strcmp(dtstr,'double') > 0 ) 
    % read data
    V = fread(fid,3*prod(sz),'double');
    V = reshape(V,[3 sz]);
    V = permute(V,[2 3 4 1]);
  end
end

fclose(fid);
    

