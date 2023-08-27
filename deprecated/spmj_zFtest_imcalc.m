function Z = spmj_zFtest_imcalc(X,c,r); 
% F-test on each voxel - one factor anova 
% function z = spmj_zFtest_imcalc(X,c,r); 
% This calculates and F-test based on the first-level betas 
if (nargin>2 && ~isempty(r)) % Runs are given, subtract the mean 
    rn=unique(r)';
    for ri=rn
        i=find(r==rn(ri)); 
        X(i,:)=bsxfun(@minus,X(i,:),sum(X(i,:))/length(i));
    end; 
    dfx=length(unique(rn));
else 
    dfx=1; 
end;

[N,P]=size(X);     % size of training set 
classes=unique(c)'; %classses we do classification on
cc=size(classes,2); %class count


% compute (estimated) mean and covariance matrix
muK=zeros([cc P]);    %means
Sw=zeros(1,P);      % Within class variability 
mu=mean(X,1);  % Overall mean 

%-------------calculate Parameter-----------------
for i=1:cc;
    j = find(c==i);                                     % select datapoints in this class
    n(i) = length(j);                                           % number of sampels per category 
    muK(i,:) = sum(X(j,:),1)/n(i);                         % get the Cluster means 
    res = bsxfun(@minus,X(j,:),muK(i,:));
    Sw = Sw+sum(res.^2,1);                         % Estimate common covariance matrix
end;
%-------------Calculate F-test----------------------
df1=cc-1;
df2=N-cc-dfx; 
Sb=bsxfun(@minus,muK,mu).^2;
Sb=bsxfun(@times,Sb,n'); 
Sb=sum(Sb,1); 
F=Sb./df1./(Sw./df2);
p=fcdf(F,df1,df2);
eps=0.000001; 
p(p>1-eps)=1-eps;
p(p<eps)=eps; 
Z=norminv(p); 
