function [D,start_sess]=spmj_get_ons_struct(SPM);
% Find the onsets for all events 
% function [D,start_sess]=spmj_get_ons_struct(SPM);
% Onsets are in image number (corrected for the 0-start in SPM!) 
% Joern Diedrichsen j.diedrichsen@ucl.ac.uk
D.block=[];         % Session number 
D.event=[];         % Event number 
D.eventname={};     % Event name 
D.num=[];           % number of event occurnace 
D.ons=[];           % Event onset in terms of the whole time series 
num_sess=length(SPM.Sess);
start_sess=1;

switch (SPM.xBF.UNITS)
    case 'scans'
        fc=1;
    case 'secs'
        fc=SPM.xY.RT;
    otherwise 
        error('unknown  units'); 
end; 

for s=1:num_sess
    % find that eventtype in that Session
    num_events=size(SPM.Sess(s).U,2);
    for e=1:num_events
            if (size(SPM.Sess(s).U(e).ons,1)>size(SPM.Sess(s).U(e).ons,2))
                event_ons=[SPM.Sess(s).U(e).ons./fc'+start_sess(s)];
            else 
                event_ons=[SPM.Sess(s).U(e).ons./fc+start_sess(s)];
            end;
            for i=1:length(event_ons)
                D.block(end+1,1)=s; 
                D.event(end+1,1)=e;
                D.eventname{end+1,1}=SPM.Sess(s).U(e).name{1};
                D.num(end+1,1)=i;
                D.ons(end+1,1)=event_ons(i);
            end;
    end;
    start_sess(s+1)=start_sess(s)+SPM.nscan(s);
end;
start_sess=start_sess(1:end-1);
