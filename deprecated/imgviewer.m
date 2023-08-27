% Copyright (C) 2000 Jarrod Millman
% Matlab verion by Joern Diedrichsen 
%
% This file is for viewing *.img brain files.
%
%
%  imgviewer (name, sx, sy, sz, slice) displays a matrix as a scaled bw image.
%
% Author: Jarrod Millman <millman@socrates.berkeley.edu>
% Created: May 2000

% This is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, write to the Free
% Software Foundation, 59 Temple Place - Suite 330, Boston, MA
% 02111-1307, USA.
function error = imgviewer (name, sx, sy, sz, slice)

  sz=[sx*sy*sz];
  slice_size=sx*sy;

  fid=fopen(name,'r');
  brain=fread(fid,sz,'int16');
  fclose(fid);
  colormap(gray);
  brighten(0.3);
  imagesc(reshape(brain(1+slice*(slice_size):(slice+1)*slice_size),sx,sy)');

  error=0;

