function [Trans,Rotat]=spmj_moveparams(SPM,varargin)
% function [Trans,Rotat]=spmj_moveparams(SPM,varargin)
% Loads and plots realignment parameters
% finds all the rp_*.txt files
% Loads first one and plots it.
% needs SPM-structure for session structure
% VARARGIN:
%   'event_rel',eventname
%   'subset',[frames]
%   'average',indicator
%   'event_lines',eventname
moveparams={};
absolute=1;
subset=1:N;
event_lines=0;
event='';
c=1;
type='plain';

while c<=length(varargin)
    switch (varargin{c})
        case 'subset'
            
        case 'abs'
            absolute=1;
            c=c+1;
        case 'dir'
            absolute=0;
            c=c+1;
        case 'event_lines'
            event=varargin{c+1};
            event_lines=1;
            c=c+2;
        case 'event_rel'
            event=varargin{c+1};
            type='event_rel';
            c=c+2;
        case 'average'
            type='average';
            indicator=varargin{c+1};
            c=c+2;
        case {'moveparams','subset'}
            assignin('caller',options{c},options{c+1});
            eval([varargin{c} '=' varargin{c+1}]);
            c=c+2;
            
        otherwise
            error(['Unknow option:' varargin{c}]);
    end;
end;

if (isempty(moveparams))
    spm_get
    P=get_files('rp_*.txt');
    Params=dlmread(P{1});
    N=size(Params,1);
    
    
    
    Trans_diff=[zeros(1,3);diff(Params(subset,1:3))];
    Trans_abs=sqrt(sum(Trans_diff.^2,2));
    [e,s]=spm_all_ons(SPM,event);
    
    Rotat_diff=[zeros(1,3);diff(Params(subset,4:6))];
    for n=1:length(subset)
        P=[zeros(1,3) Rotat_diff(n,:) ones(1,3) zeros(1,3)];
        M=spm_matrix(P);
        b=M*[1;0;0;1];
        Rotat_abs(n,:)=acos([1;0;0]'*b(1:3));
    end;
    
    f=180/pi;
    t=1:length(subset);
    Trans=Params(subset,1:3);
    Rotat=Params(subset,4:6)*f;
    
    
    switch (type)
        case 'plain'
            subplot(2,1,1);
            ts_plot(t,Trans,'Translation [mm]',s,event_lines,e,{'x','y','z'});
            subplot(2,1,2);
            ts_plot(t,Rotat,'Rotation [deg]',s,event_lines,e,{'Pitch','Roll','Yaw'});
            
        case 'event_rel'
            subplot(4,1,1);
            ts_plot(t,Trans,'Translation [mm]',s,event_lines,e,{'x','y','z'});
            subplot(4,1,2);
            ts_plot(t,Rotat,'Rotation [deg]',s,event_lines,e,{'Pitch','Roll','Yaw'});
            t_ev=[-1:7];
            [TER,sd_TER]=evoked_response(Trans_abs,e,t_ev);
            [RER,sd_RER]=evoked_response(Rotat_abs,e,t_ev);
            subplot(4,1,3);
            plot(t_ev,TER,'r-',t_ev,TER+sd_TER,'r:',t_ev,TER-sd_TER,'r:');
            ylabel('Translation');
            xlabel('Image');
            subplot(4,1,4);
            plot(t_ev,RER*f,'r-',t_ev,(RER+sd_RER)*f,'r:',t_ev,(RER-sd_RER)*f,'r:');
            ylabel('Rotation [deg]');
            xlabel('Image');
        case 'average'
            subplot(4,2,[1:2]);
            ts_plot(t,Trans,'Translation [mm]',s,event_lines,e,{'x','y','z'});
            subplot(4,2,[3:4]);
            ts_plot(t,Rotat,'Rotation [deg]',s,event_lines,e,{'Pitch','Roll','Yaw'});
            CAT.markercolor={'b','g','r'};
            CAT.markertype={'x','v','*'};
            subplot(4,2,[5 7]);
            Trans=Trans-repmat(mean(Trans),size(Trans,1),1);
            Rotat=Rotat-repmat(mean(Rotat),size(Rotat,1),1);
            [TM,TSE]=lineplot(Trans,indicator,'CAT',CAT,'markersize',8,'legend',{'x','y','z'});
            ylabel('Translation (mm)');
            subplot(4,2,[6 8]);
            [RM,RSE]=lineplot(Rotat,indicator,'CAT',CAT,'markersize',8,'legend',{'Pitch','Roll','Yaw'});
            ylabel('Rotation (deg)');
            
    end;
end;
function ts_plot(t,Y,yl,s,event_lines,e,leg);
plot(t,Y);
ylabel(yl);
xlabel('Image');
drawlines(s,[0 0 0]);
if (event_lines)
    drawlines(e,[1 0 0]);
end;
legend(leg);