function error = spmj_flipimage(filename,varargin)
% [A,hdr]=flipimage (name, varargin)
%   uses spm_read_hdr to asess the header information
%   filename can be either a single filename
%   or a filename including a placeholder 'name*.img'
%   varargin can be
%   'ud' flips the image in y-direction
%   'lr' flips the image in x-direction
%   'tb' flips the image in z-direction
c=1;
isUD=0;scY=1;
isLR=0;scX=1;
isTB=0;scZ=1;
outfilename=[]; 
while c<=length(varargin)
    switch (varargin{c})
        case {'ud','y'}
            isUD=1;scY=-1;c=c+1;
        case {'lr','x'}
            isLR=1;scX=-1;c=c+1;
        case {'tb','z'}
            isTB=1;scZ=-1;c=c+1;
        case 'outfilename'
            outfilename=varargin{c+1}; 
            c=c+2;
    end;
end;
[dir,name,ext]=spm_fileparts(filename);
if (isempty(outfilename))
    outfilename=fullfile(dir,[name ext]);
end;

VS=spm_vol(filename);
for n=1:length(VS)
    VT(n,1)    = struct('fname',	outfilename,...
        'dim',		VS(n).dim,...
        'dt',       VS(n).dt,...
        'mat',		VS(n).mat*spm_matrix([isLR*(VS(n).dim(1)+1) isUD*(VS(n).dim(2)+1) isTB*(VS(n).dim(3)+1) 0 0 0 scX scY scZ]),...
        'pinfo',	VS(n).pinfo,...
        'n',        VS(n).n,...
        'descrip',	'flipped');
end;

VT    = spm_create_vol(VT);
X=spm_read_vols(VS);
if (isUD==1)
    X=flipdim(X,2);
end;
if (isLR==1)
    X=flipdim(X,1);
end;
if (isTB==1)
    X=flipdim(X,3);
end;

for n=1:length(VT) 
    for z=1:VT(n).dim(3)
        VT(n)=spm_write_plane(VT(n),X(:,:,z,n),z);  % Avoid scaling of volume 
    end;
end;


