function varargout=spmj_crossval_std(xX,Y,Sess); 
% Standard estimation and crossvalidation: 
% Plugin function into spmj_spm_crossval
if (nargin==0)
    varargout={2,1};         % Number of betaStats images 
    return; 
end; 
if (isempty(Y))      % Do precomputations on the Noise regressor only matrix for each run seperately
    for i=1:max(xX.iiS)
        xN=xX.xKXs.X(Sess(i).row,xX.iiN>0 & xX.iiS==i);
        xX.xXNs(i)=spm_sp('Set',xN); 
        xX.pXN{i}=spm_sp('x-',xX.xXNs(i)); 
    end; 
    varargout={xX}; 
    return; 
end; 


beta  = xX.pKX*Y;                    %-Parameter estimates
res   = spm_sp('r',xX.xKXs,Y);       %-Residuals
ResSS(1,:) = sum(res.^2);                   %-Residual SSQ

nSess=length(Sess);
RSS=zeros(3,size(Y,2));
BRSS=zeros(2,size(Y,2)); 
BTSS=zeros(2,size(Y,2)); 
beta_hat=zeros(size(beta));
for i=1:nSess
    % Get the prediction from the 
    y=Y(Sess(i).row,:);   % Get data from this block 
    j=find(xX.iiS==i & xX.iiC>0);  % Find betas relevant for this block 
    RSS(1,:)=RSS(1,:)+sum(bsxfun(@minus,y,mean(y)).^2); 
    if (any(xX.iiN(xX.iiS==i)==1))          % Any nuisance regressors? 
        res=spm_sp('res',xX.xXNs(i),y);
        RSS(2,:)=RSS(2,:)+sum(res.^2); 
    else 
        RSS(2,:)=RSS(1,:);
    end; 
    
    % Make the transfer betas from the other blocks 
    for k=1:length(j) 
        beta_hat(j(k),:)=mean(beta(xX.iiS~=i & xX.iiC==xX.iiC(j(k)) & xX.iiT==xX.iiT(j(k)) & xX.iiB==xX.iiB(j(k)),:));         
    end; 
    
    yhat=xX.xKXs.X(Sess(i).row,j)*beta_hat(j,:); 
    yres=y-yhat; 
    res=spm_sp('res',xX.xXNs(i),yres);
    RSS(3,:)=RSS(3,:)+sum(res.^2); 
    
    for c=1:2 
        j=find(xX.iiS==i & xX.iiC==c); 
        mB=bsxfun(@minus,beta(j,:),mean(beta(j,:))); 
        mBh=bsxfun(@minus,beta_hat(j,:),mean(beta_hat(j,:))); 
        BTSS(c,:)=BTSS(c,:)+sum(mB.^2);
        BRSS(c,:)=BRSS(c,:)+sum((mB-mBh).^2);
    end;
end; 
betaStats=(BTSS-BRSS)./BTSS;   % Beta stats is the amount of variance predicted of the betas 
varargout={beta,RSS,betaStats,ResSS};

