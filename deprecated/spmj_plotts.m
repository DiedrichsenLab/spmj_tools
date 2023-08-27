function [t,y]=spmj_plotts(varargin)
% Needs the following variables in workspace: 
%   SPM: main SPM structure  
% Which voxels: 
%   'ROI',filename: 
%       file describing a ROI, over which the function is plotted
%   'xyz': 
%       x,y,z cordindiantes to extract the values from 
% Type of plot 
%   'event_lines'
%   'event_rel',name,num_frames:   display an average event related response+
%                       the fitted version from the basis functions 
% Other 
%   'figure',['spm'/'seperate']
% Defaults 
Col   = [0 0 0; .8 .8 .8; 1 .5 .5];
event='';
type='plain';
data='xyz';
fig='seperate';
% Deal with user input 
% ----------------------------------------------------------------------
c=1;
while c<=length(varargin)
    switch varargin{c}
        case 'ROI'
            data='ROI';
            roi_name=varargin{c+1};
            c=c+2; 
        case 'scan'
            scan=varargin{c+1}
            c=c+2;
        case 'event_lines'
            event=varargin{c+1};
            c=c+2;
        case {'event_rel'}
            type='event_rel';
            event=varargin{c+1};
            num_frames=varargin{c+2};
            c=c+3;
        case 'figure'
            fig=varargin{c+1};
            c=c+2;
        otherwise 
            error(['Unknown option: ' varargin{c}]);
    end;
end;


try 
    SPM=evalin('caller','SPM');
catch
    error('Structure SPM was not given in calling workspace');
end;

if (data=='ROI')
    if (isempty(roi_name))
        roi_name=spm_get(1,'*.img','Select ROI file');
    end;
    X=spm_read_vols(spm_vol(roi_name));
    [x,y,z]=ind2sub(size(X),find(X>0));
    XYZ=[x';y';z';ones(1,length(x))];
    y_raw = spm_get_data(SPM.xY.VY,XYZ);
    y_raw = mean(y_raw,2);
    y = spm_filter(SPM.xX.K,SPM.xX.W*y_raw);
    XYZstr = ['ROI:' roi_name];    
end;

% recompute the model, residuals, BCOV and adjusted data
%  ---------------------------------------------------------------------
beta=SPM.xX.pKX*y;   % pKX = inv ((WX)'*(WX))*(WX)'*Wy
R=y-SPM.xX.xKXs.X*beta;
ResMS=R'*R/SPM.xX.trRV;
Bcov  = ResMS*SPM.xX.Bcov;
indexADJ=[SPM.xX.iB SPM.xX.iG];
indexInterest=SPM.xX.iC;
y_adj=y-SPM.xX.xKXs.X(:,indexADJ)*beta(indexADJ);
y_pred=SPM.xX.xKXs.X(:,indexInterest)*beta(indexInterest);
CI=1.65;

% Extract Session Structure 
% -------------------------------------------------------
num_sess=length(SPM.Sess);
it=[1:length(y_adj)];is=[];isess=[];
for s=1:num_sess
    is=[is 1:SPM.nscan(s)];
    isess=[isess ones(1,SPM.nscan(s))*s];
end;
TITLE='Time series';

[event_ons,start_sess]=spmj_all_ons(SPM,event);

% --------------------------------------------------------------
%-Get Graphics figure handle
current_fig=gcf;
Fgraph = spm_figure('GetWin','Graphics');
switch (fig)
    case 'seperate'
        if (mod(current_fig,1)~=0 | current_fig==Fgraph)
            figure;
        else 
            figure(current_fig);
            cla;
        end;
    case 'spm'
        %-Delete previous axis and their pagination controls (if any)
        spm_results_ui('Clear',Fgraph,2);
        figure(Fgraph);
        subplot(2,1,2);
        cla;
end;
hold on;

switch type 
    case 'plain'
        plot(it,y_adj,'Color',[0 0 0]);
		plot(it,y_pred,'LineWidth',4,'Color',[0.4 0.4 0.4]);

        % Scan borders
        drawlines(it(find(is==1)),'k');

        % Events in lines 
        if (~isempty(event_ons))
            drawlines(event_ons,'k');
        end;

        title(TITLE,'FontSize',12)
        xlabel('Scan')
        ylabel(['BOLD: ',XYZstr])
    case 'event_rel'
        [y_ev,y_sd]=evoked_response(y_adj,event_ons,num_frames);
        [y_pev,y_psd]=evoked_response(y_pred,event_ons,num_frames);
        % calculate CI and fitted evoked response in high resolution
        dt=SPM.xBF.dt;
        X     = SPM.xBF.bf/dt;
        C=SPM.xCon(1).c';
        B=C*beta;     
        y_fit=X*B;
        y_fit_ci = CI*sqrt(diag(X*C*Bcov*C'*X'));

        time_low   = [0:num_frames-1]*SPM.xY.RT;
        % time_high  = ([1:size(X,1)] - 1)*SPM.xBF.dt;
        T=length(event_ons); % number of events 
        errorbar(time_low,y_ev,y_sd/sqrt(T)*CI);
		plot(time_low,y_ev,'LineWidth',4,'Color',[0.4 0.4 0.4]);
		% plot(time_high,y_fit,'r-',time_high,y_fit+y_fit_ci,'r.',time_high,y_fit-y_fit_ci,'r.');
		plot(time_low,y_pev,'r-',time_low,y_pev+y_psd,'r.',time_low,y_pev-y_psd,'r.');
        xlabel('Time(s)')
        ylabel(['BOLD: ',XYZstr]);
        % resample at high frquency 
        %t=time_high;
        % y=[y_fit interp1(time_low,y_ev,time_high,'pchip')'];
end;        
hold off;
%-call Plot UI
%----------------------------------------------------------------------
%spm_results_ui('PlotUi',gca)
