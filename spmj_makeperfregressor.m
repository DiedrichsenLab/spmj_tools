function X=spmj_makeperfregressor(N,onset,dur,varargin); 
% Make a perfusion regression model based on 
% 1.Simple BOLD regressor 
% 2. Perfusion baseline regressor 
% 3. Functional Perfusion regressor Perf x BOLD 
TR=2.0; 
bf=spm_hrf(2); 
X=zeros(N,3); 
if (length(onset)~=length(dur))
    error ('duration and onset need to be same length'); 
end; 
for i=1:length(onset) 
    X(onset(i):onset(i)+dur(i))=1;
end; 
% plot(X(:,1)); hold on; 
X(:,1)=spmj_conv(X(:,1),bf); 
% plot(X(:,1)); 
X(:,2)=repmat([1;-1],N/2,1); 
X(:,3)=X(:,1).*X(:,2); 
