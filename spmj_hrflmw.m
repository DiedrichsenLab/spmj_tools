function y=spmj_hrflmw(t,P,hrf)
% function y=spmj_hrflmw(t,P,hrf)
% Hemodynamic response function with parameters PS
% convolved with a boxcar function of magnitude,latency and width (1 s)
dt=1/8;
T=[-3:dt:32]';
B=zeros(length(T),1);
[dummy,i1]=min(abs(P(1)-T));
[dummy,i2]=min(abs(P(3)+P(1)-T));
B(i1:i2)=P(2);
s=((T(i1)-P(1))+dt/2)/dt;
if (s<0 | s>1) 
    error('AR');
end;
B(i1)=s*P(2);
s=((T(i2)-P(1)-P(3))+dt/2)/dt;
if (s<0 | s>1) 
    error('AR');
end;
B(i2)=(1-s)*P(2);
H=conv(B,hrf);
H=H(1:length(T));
y=interp1(T,H,t);