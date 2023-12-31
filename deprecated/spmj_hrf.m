function [hrf,p] = spm_hrf(RT,P,T)
% Return a hemodynamic response function
% FORMAT [hrf,p] = spm_hrf(RT,[p])
% RT   - scan repeat time
% p    - parameters of the response function (two gamma functions)
%
%                                                     defaults
%                                                    (seconds)
%   p(1) - delay of response (relative to onset)         6
%   p(2) - delay of undershoot (relative to onset)      16
%   p(3) - dispersion of response                        1
%   p(4) - dispersion of undershoot                      1
%   p(5) - ratio of response to undershoot               6
%   p(6) - onset (seconds)                               0
%   p(7) - duration of underlying neural event (s)       0
%   p(8) - length of kernel (seconds)                   32
%
% T    - microtime resolution [Default: 16]
%
% hrf  - hemodynamic response function
% p    - parameters of the response function
%__________________________________________________________________________
% Copyright (C) 1996-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hrf.m 4472 2011-09-08 17:42:32Z guillaume $


%-Default parameters
%--------------------------------------------------------------------------
p   = [6 16 1 1 6 0 0 32];
if nargin > 1
    p(1:length(P)) = P;
end

%-Microtime resolution
%--------------------------------------------------------------------------
if nargin > 2
    fMRI_T = T;
else
    fMRI_T = spm_get_defaults('stats.fmri.t');
end

%-Modelled hemodynamic response function - {mixture of Gammas}
%--------------------------------------------------------------------------
dt  = RT/fMRI_T;
u   = [0:(p(8)/dt)] - p(6)/dt;
hrf = spm_Gpdf(u,p(1)/p(3),dt/p(3)) - spm_Gpdf(u,p(2)/p(4),dt/p(4))/p(5);
neur = zeros(size(u)); 
indx=find((u*dt)<p(7)); 
neur(indx)=1;
nextIndex=indx(end)+1; % Index after 
neur(nextIndex)=(u(nextIndex)*dt-p(7))/dt; 
hrf=conv(hrf,neur);
hrf = hrf([0:(p(8)/RT)]*fMRI_T + 1);
hrf = hrf'/sum(hrf);
