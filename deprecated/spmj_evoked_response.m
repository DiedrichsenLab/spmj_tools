function [m,sd,ERP]=evoked_response(Y,tau,leng)
% function [m,sd,ERP]=evoked_response(Y,tau,leng,fig)
% extracts and displays evoked response by simple averaging 
% INPUTS: 
%   Y: a column per p time series
%   tau: time of events 
%   leng: length of the extracted time series
%         scalar: 0-length-1
%         vector: indicates the frames, e.g. [-2:7]
% OUTPUTS:
%   m: mean evoked response
%   sd: standard deviation of evoked response 
%   ERP: cell array of individual responses
[N,P]=size(Y);m=[];ERP={};
if (length(leng)==1)
    frames=[0:leng-1];
else
    frames=leng;
    leng=length(frames);
end;
for p=1:P 
    for t=1:length(tau)
        r=ones(leng,1)*NaN;
        r(1:leng,1)=save_index(Y(:,p),tau(t)+frames);
        ERP{p}(t,:)=r;
    end;
    m(p,:)=nanmean(ERP{p});
    sd(p,:)=nanstd(ERP{p});
end;

