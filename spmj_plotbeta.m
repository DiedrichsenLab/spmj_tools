function spm_plotbeta(Ic,varargin)
% Needs the following variables in workspace: 
%   SPM: main SPM structure  xSPM: specific selected contrast hReg: handle
%   to cursor pointer 
% Optional: 
%   xY: structure describing a ROI, over which the function is plotted

%-Get Graphics figure handle
Fgraph = spm_figure('GetWin','Graphics');
%-Delete previous axis and their pagination controls (if any)
spm_results_ui('Clear',Fgraph,2);

Col   = [0 0 0; .8 .8 .8; 1 .5 .5];

% - try to get SPM, xSPM and hReg 
% ----------------------------------------------------------------------
% prevent 'SPM' from restarting SPM
%if (eval('caller','exist(''SPM'')')==2)
%    error('Structure SPM is not given in calling workspace');
%end;    
try 
    SPM=evalin('caller','SPM');
    xSPM=evalin('caller','xSPM');
    hReg=evalin('caller','hReg');
catch
    error('Structure xSPM, or hReg were not given in calling workspace');
end;
try 
    VOI=evalin('caller','xY');
catch
    VOI=[];
end;

% ----------------------------------------------------------------------
% Extract relevant statistics and times-series from the voxel or VOI
% 
% ----------------------------------------------------------------------
%-Find nearest voxel [Euclidean distance] in point list & update GUI
%-----------------------------------------------------------------------
if (isempty(VOI))
    if ~length(xSPM.XYZmm)
        spm('alert!','No suprathreshold voxels!',mfilename,0);
        Y = []; y = []; beta = []; Bcov = [];
        return
    end
    [xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
    spm_XYZreg('SetCoords',xyz,hReg);
    XYZ     = xSPM.XYZ(:,i);		% coordinates
    
    %-Extract filtered and whitened data from files
    %======================================================================
    %=
    try
        y_raw = spm_get_data(SPM.xY.VY,XYZ);
        y = spm_filter(SPM.xX.K,SPM.xX.W*y_raw);
    catch
        try
            % remap files in SPM.xY.P if SPM.xY.VY is no longer valid
            %-------------------------------------------------------
            SPM.xY.VY = spm_vol(SPM.xY.P);
            y_raw = spm_get_data(SPM.xY.VY,XYZ);
            y = spm_filter(SPM.xX.K,SPM.xX.W*y_raw);
        catch
            error('Original data have been moved or renamed');
        end
    end
    XYZstr = sprintf(' at [%g, %g, %g]',xyz);
else
    XYZ=SPM.xVol.iM*[VOI.XYZmm;ones(1,size(VOI.XYZmm,2))];
    y_raw = spm_get_data(SPM.xY.VY,XYZ);
    y_raw = mean(y_raw,2);
    y = spm_filter(SPM.xX.K,SPM.xX.W*y_raw);
    XYZstr = ['VOI:' VOI.name];    
end;

% recompute the model, residuals, BCOV and adjusted data
%  ---------------------------------------------------------------------
beta=SPM.xX.pKX*y;
R=y-SPM.xX.xKXs.X*beta;
ResMS=R'*R/SPM.xX.trRV;
Bcov  = ResMS*SPM.xX.Bcov;
indexADJ=[SPM.xX.iB SPM.xX.iG];
y_adj=y-SPM.xX.xKXs.X(:,indexADJ)*beta(indexADJ);
num_sess=length(SPM.Sess);
CI=1.65;

% Get the contrasts of interest 
if (nargin==0 | isempty(Ic))
    Ic    = spm_input('Which contrast?','!+1','m',{SPM.xCon.name});
end;

TITLE = {SPM.xCon(Ic).name};
if xSPM.STAT == 'P'
    TITLE = {Cplot SPM.xCon(Ic).name '(conditional estimates)'};
end;

% Deal with user input 
% ----------------------------------------------------------------------
c=1;
while c<=length(varargin)
    switch varargin{c}
        case {'Contrast'}
    end;
end;


% compute contrast of parameter estimates and 90% C.I.
cbeta = SPM.xCon(Ic).c'*beta;
CI    = 1.6449;					% = spm_invNcdf(1 - 0.05);
CI    = CI*sqrt(diag(SPM.xCon(Ic).c'*Bcov*SPM.xCon(Ic).c));
% bar chart
figure(Fgraph)
subplot(2,1,2)
cla
hold on

% estimates
h     = bar(cbeta);
set(h,'FaceColor',Col(2,:))

% standard error
%--------------------------------------------------------------
for j = 1:length(cbeta)
    line([j j],([CI(j) 0 - CI(j)] + cbeta(j)),...
        'LineWidth',6,'Color',Col(3,:))
end

title(TITLE,'FontSize',12)
xlabel('contrast')
ylabel(['contrast estimate',XYZstr])
set(gca,'XLim',[0.4 (length(cbeta) + 0.6)])
hold off

%-call Plot UI
%----------------------------------------------------------------------
spm_results_ui('PlotUi',gca)
