function spmj_plot_Vi(SPM,varargin); 
% Series of diagnostic graphs for the estimated covariance structure of an SPM.  
% and for the consequential prewhitening. Also takes into account the
% filtering. Current version assumes that all runs have the same lenght. 
numRuns = length(SPM.nscan); 
numImages = SPM.nscan(1); 
if std(SPM.nscan)>0 
    error('all runs have to have the same length'); 
end; 

for r=1:numRuns
    indx = SPM.Sess(r).row; 
    V=full(SPM.xVi.V(indx,indx)); 
    W=full(SPM.xX.W(indx,indx)); 
    Vth=pinv(W*W'); 
    KX=SPM.xX.K(1).X0; 
    K=(eye(numImages)-KX*pinv(KX));
    KVK=K*V*K'; 
    subplot(3,2,1); 
    imagesc(pinv(K)*V*pinv(K)); 
    subplot(3,2,2); 
    imagesc(pinv(K*W*W'*K')); 
    
    
    
    
    
    tsV(r,:)=V(1,:); 
    tsVth(r,:)=Vth(1,:); 
    tsW(r,:)=W(1,:); 
    tsK(r,:)=K(1,:); 
    tsKVK(r,:)=KVK(1,:); 
end; 
xlim=[0 size(tsV,2)]; 
t=[1:numImages]; 
subplot(3,2,1); 
traceplot(t,tsV,'errorfcn','stderr'); 
set(gca,'xlim',xlim); 
subplot(3,2,2); 
plot(t,tsVth); 
set(gca,'xlim',xlim); 
subplot(3,2,3); 
plot(t,tsW); 
set(gca,'xlim',xlim); 
subplot(3,2,4); 
plot(t,tsK); 
set(gca,'xlim',xlim); 



