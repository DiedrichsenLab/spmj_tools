function [W,TOA]=exponential_iti(meanTOA,minTOA,N)
% provides an approximation for exponentially distributed Inter-trial-intervals 
% with a minimum ITT of  x
% for a run of N
% then it takes each integer-step
x=[minTOA-1:N]';
p=geopdf(x,1/meanTOA);
num_events=N/meanTOA;
W=round(p.*num_events/sum(p));
TOA=x+1;
indx=find(W==0);
W(indx)=[];
TOA(indx)=[];
% now equate that the sum of TOA's is exactly N
while  (N-W'*TOA~=0)
    x=N-W'*TOA;
    i=find(TOA==abs(x));
    if (isempty(i))
         i=unidrnd(length(TOA));
    end;
    W(i)=W(i)+1*sign(x);
    indx=find(W==0);
    W(indx)=[];
    TOA(indx)=[];
end;
