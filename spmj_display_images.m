function [minval,maxval]=spm_display_images(P,varargin)
% function spm_display_images(Images,varargin)
%   options:  
%       slices: which slices to present 
%       orient: default:  
%       rows: number of rows
%       columns: number of columns
%       

slices=[];
orient='columns'; 
colorbar='none';
dim=[];
c=1; 
colorm='gray';
minval=[];maxval=[];
while c<length(varargin)
    switch (varargin{c})
        case {'slices','orient','colorbar','minval','maxval'}
            eval([ varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option :%s',varargin{c}));
    end;
end;
try
    Fgraph= spm_figure('FindWin','Graphics');
    spm_results_ui('Clear',Fgraph,2);
    figure(Fgraph);
catch
    figure(4);
end;

if (isempty(P)) 
    P=spm_get(inf,'*.img','Images to present');
end;
if (iscell(P))
    P=char(P);
end;

% calculate rows/columns/slices
if (isempty(slices))
    slices=[1:V(1).dim(3)];
end;
if (strcmp(orient,'columns'))
    rows=length(slices);
    columns=size(P,1);
else
    rows=size(P,1);
    columns=length(slices);
end;
    
V=spm_vol(P);
if (isempty(dim))
    dim=V(1).dim(1:2);
end;
for i=1:size(P,1)
    for j=1:length(slices);
        M=eye(4);
        M(3,4)=slices(j);
        X(:,:,j,i)=spm_slice_vol(V(i),M,dim,1);
    end;
end;
if (isempty(maxval));
    maxval=max(X(:));
end;
if (isempty(minval))
    minval=min(X(:));
end;
for i=1:size(X,3)
    for j=1:size(X,4)
        subplot(rows,columns,(i-1)*columns+j);
        imagesc(flipud(X(:,:,i,j)'),[minval maxval]);
        axis square
        if (i==1) 
            dash=findstr(P(j,:),'\');
            title(P(j,dash(end)+1:end));
        end;
        if(j==1)
            ylabel(sprintf('Slice %d',slices(i)));    
        end;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
    end;
end;
colormap(colorm);
