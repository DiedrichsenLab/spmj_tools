function varargout=spmj_crossval_ridge(xX,Y,Sess); 
% estimation considering all regression coefficients as random effects:
% Estimation of the hyper parameters happens for each scan independently 
% Plugin function into spmj_spm_crossval
if (nargin==0)
    varargout={2,1};         % Number of betaStats images and of hyperparamters (msRes)
    return; 
end; 
if (isempty(Y))      % Do precomputations on the Noise regressor only matrix 
    numCond=max(xX.iiC); 
    xN=xX.xKXs.X(:,xX.iiN>0);
    xX.xXNs=spm_sp('Set',xN); 
    xX.pXN=spm_sp('x-',xX.xXNs); 
    % Now do the regressors of interest 
    xI=xX.xKXs.X(:,xX.iiN==0);
    xX.xI.iiC=xX.iiC(xX.iiN==0); 
    xX.xI.iiT=xX.iiT(xX.iiN==0); 
    xX.xI.iiB=xX.iiB(xX.iiN==0); 
    xX.xI.iiS=xX.iiS(xX.iiN==0); 
    
    xX.xI.X=spm_sp('res',xX.xXNs,xI); 
    xX.xXIs=spm_sp('Set',xX.xI.X); 
    xX.pXI=spm_sp('x-',xX.xXIs); 
    varargout={xX}; 
    return;
end; 


alphaVal=[0 exp([-4:0.5:4])]; 
numAlpha=length(alphaVal); 

nSess=length(Sess);
RSS=zeros(3,size(Y,2),numAlpha);
BRSS=zeros(2,size(Y,2),numAlpha); 
BTSS=zeros(2,size(Y,2),numAlpha); 

% Ridge parameter is only applied to regressors of interest 
beta=zeros(size(xX.xXIs.X,2),size(Y,2),numAlpha)*NaN; 
beta_hat=zeros(size(beta));

% get all betas for different values of the ridge parameter 
sX=xX.xXIs; 
r=size(sX.ds,1); 
for a=1:numAlpha
    beta(:,:,a)=sX.v(:,1:r)*diag( sX.ds(1:r)./(sX.ds(1:r).^2+alphaVal(a).^2))*sX.u(:,1:r)'*Y; 
end;


% First, establish the RSS after subtract of all noise regressors
r=spm_sp('res',xX.xXNs,Y); 
RSS(2,:,1)=sum(r.^r); 


% Check crossvalidation 
for i=1:nSess

    j=find(xX.xI.iiS==i); 
    nj=find(xX.xI.iiS~=i); 

    % Make the transfer betas from the other blocks 
    
    for a=1:numAlpha
        for l=1:length(j) 
            beta_hat(j(l),:,a)=mean(beta(xX.iiC==xX.iiC(j(l)) & xX.iiT==xX.iiT(j(l)) & xX.iiB==xX.iiB(j(l)),:,a));         
        end; 
        y=r(Sess(i).row,:);   % Get data from this block 
        yhat=xX.xI.X(Sess(i).row,j)*beta_hat(j,:,a); 
        yres=y-yhat; 
        RSS(3,:,a)=RSS(3,:,a)+sum(yres.^2); 
    end; 
    
    % Check regressor RSS 
    for c=1:2 
        k=find(xX.xI.iiS==i & xX.xI.iiC==c); 
        for a=1:numAlpha 
            mB=bsxfun(@minus,beta(k,:,a),mean(beta(k,:,a))); 
            mBh=bsxfun(@minus,beta_hat(k,:,a),mean(beta_hat(k,:,a))); 
            BTSS(c,:,a)=BTSS(c,:,a)+sum(mB.^2);
            BRSS(c,:,a)=BRSS(c,:,a)+sum((mB-mBh).^2);
        end; 
    end;
end; 
betaStats=(BTSS-BRSS)./BTSS;   % Beta stats is the amount of variance predicted of the betas 
b2=squeeze(betaStats(2,:,:));
plot(b2');
varargout={beta,RSS,betaStats,mean(h,3)}; 