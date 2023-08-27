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
RSS=zeros(3,size(Y,2));
BRSS=zeros(2,size(Y,2));
BTSS=zeros(2,size(Y,2));

% Ridge parameter is only applied to regressors of interest
beta=zeros(size(xX.X,2),size(Y,2))*NaN;

betaA=zeros(size(xX.xXIs.X,2),size(Y,2),numAlpha)*NaN;
beta_hat=zeros(size(xX.xXIs.X,2),size(Y,2))*NaN;
beta_h=zeros(size(betaA));

% get all betas for different values of the ridge parameter
sX=xX.xXIs;
r=size(sX.ds,1);
for a=1:numAlpha
    betaA(:,:,a)=sX.v(:,1:r)*diag( sX.ds(1:r)./(sX.ds(1:r).^2+alphaVal(a).^2))*sX.u(:,1:r)'*Y;
end;

% First, establish the RSS after subtract of all noise regressors
indxNBeta=find(xX.iiN>0); 
indxIBeta=find(xX.iiN==0); 

beta(indxNBeta,:)=xX.pXN*Y; 
r=spm_sp('res',xX.xXNs,Y);
RSS(2,:)=sum(r.^2);

regtypes=unique([xX.xI.iiC' xX.xI.iiT' xX.xI.iiB'],'rows'); 

% For the regression leave one of the runs out
for i=1:nSess
    otherSessions=[1:nSess];
    otherSessions(i)=[];
    rssC=zeros(numAlpha,size(Y,2));
    for j=otherSessions
        target=(xX.xI.iiS==j);
        source=(xX.xI.iiS~=j & xX.xI.iiS~=i);
        
        y=r(Sess(j).row,:);   % Get data from this block
        
        for l=size(regtypes,1) 
            indS=find(source & xX.xI.iiC==regtypes(l,1) & xX.xI.iiT==regtypes(l,2) & xX.xI.iiB==regtypes(l,3)); 
            indT=find(target & xX.xI.iiC==regtypes(l,1) & xX.xI.iiT==regtypes(l,2) & xX.xI.iiB==regtypes(l,3)); 
            beta_h(indT,:,:)=repmat(sum(betaA(indS,:,:),1)./length(indS),[length(indT) 1 1]);
        end; 
            
        for a=1:numAlpha
            yhat=xX.xI.X(Sess(j).row,target)*beta_h(target,:,a);
            yres=y-yhat;
            rssC(a,:)=rssC(a,:)+sum(yres.^2);
        end;
    end;
    
    % Now select ridge parameter that brings us the lowest crossvalidation
    % error
    [~,minAlpha(i,:)]=min(rssC);
    target=(xX.xI.iiS==i);
    source=(xX.xI.iiS~=i);
    
    for l=1:size(regtypes,1)
        indS=find(source & xX.xI.iiC==regtypes(l,1) & xX.xI.iiT==regtypes(l,2) & xX.xI.iiB==regtypes(l,3)); 
        indT=find(target & xX.xI.iiC==regtypes(l,1) & xX.xI.iiT==regtypes(l,2) & xX.xI.iiB==regtypes(l,3)); 
        
        for a=unique(minAlpha(i,:))
            v=minAlpha(i,:)==a; 
            beta_hat(indT,v)=repmat(sum(betaA(indS,v,a),1)./length(indS),[length(indT) 1 1]);
        end;
    end;
    y=r(Sess(i).row,:);   % Get data from this block
    yhat=xX.xI.X(Sess(i).row,target)*beta_hat(target,:);
    yres=y-yhat;
    RSS(3,:)=RSS(3,:)+sum(yres.^2);
    
    % find the best guesses for beta
    for a=unique(minAlpha(i,:))
        k=find(xX.xI.iiS==i); 
        v=minAlpha(i,:)==a; 
        beta(indxIBeta(k),v)=betaA(k,v,a); 
    end; 
    
    for c=1:2
        k=find(xX.xI.iiS==i & xX.xI.iiC==c);
        mB=bsxfun(@minus,beta(indxIBeta(k),:),mean(beta(indxIBeta(k),:)));
        mBh=bsxfun(@minus,beta_hat(k,:),mean(beta_hat(k,:)));
        BTSS(c,:)=BTSS(c,:)+sum(mB.^2);
        BRSS(c,:)=BRSS(c,:)+sum((mB-mBh).^2);
    end; 
end;
betaStats=(BTSS-BRSS)./BTSS;   % Beta stats is the amount of variance predicted of the betas
varargout={beta,RSS,betaStats,mean(minAlpha)};