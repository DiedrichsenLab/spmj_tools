function varargout=spmj_crossval_rfx(xX,Y,Sess); 
% estimation considering all regression coefficients as random effects:
% Estimation of the hyper parameters happens for each scan independently 
% Plugin function into spmj_spm_crossval
if (nargin==0)
    varargout={2,3};         % Number of betaStats images and of hyperparamters (msRes)
    return; 
end; 
if (isempty(Y))      % Do precomputations on the Noise regressor only matrix for each run seperately
    numCond=max(xX.iiC); 
    for i=1:max(xX.iiS)
        xN=xX.xKXs.X(Sess(i).row,xX.iiN>0 & xX.iiS==i);
        xX.xXNs(i)=spm_sp('Set',xN); 
        xX.pXN{i}=spm_sp('x-',xX.xXNs(i)); 
        
        for c=1:numCond       % Condition covariances 
            for t=1:5 
                xC{c}(:,t)=sum(xX.xKXs.X(Sess(i).row,xX.iiN==0 & xX.iiC==c & xX.iiS==i & xX.iiT==t),2); 
            end; 
            xX.Q{i,c}=xC{c}*xC{c}';
            j=find(xX.iiS==i & xX.iiN<2);   % get the relevant matrix  
            xX.Gc{i,c}=find(xX.iiN(j)==0 & xX.iiC(j)==c); 
        end; 
        xX.Q{i,numCond+1}=xN*xN';
        xX.Gc{i,numCond+1}=find(xX.iiN(j)==1); 
        xX.Q{i,numCond+2}=speye(length(Sess(i).row)); 
    end; 
    
    % Scale Q's 
    N=size(xX.Q{1,1},1); 
    for i = 1:size(xX.Q,1);
        for j=1:size(xX.Q,2); 
            xX.sh(i,j) = trace(xX.Q{i,j})/length(Sess(i).row);
            xX.Q{i,j}  = xX.Q{i,j}/xX.sh(i,j);
            xX.Q{i,j}(abs(xX.Q{i,j})<0.0001)=0;
            a=sum(xX.Q{i,j}(:)~=0)./(N*N);
            if (a) <0.2 
                xX.Q{i,j}=sparse(xX.Q{i,j});
            else 
                xX.Q{i,j}=xX.Q{i,j};
            end; 
                
        end; 
    end
    varargout={xX}; 
    return; 
end; 


% beta  = xX.pKX*Y;                    %-Parameter estimates
% res   = spm_sp('r',xX.xKXs,Y);       %-Residuals
% ResSS(1,:) = sum(res.^2);                   %-Residual SSQ

nSess=length(Sess);
RSS=zeros(3,size(Y,2));
BRSS=zeros(2,size(Y,2)); 
BTSS=zeros(2,size(Y,2)); 
beta=zeros(size(xX.X,2),size(Y,2))*NaN; 
beta_hat=zeros(size(beta));
for i=1:nSess
    j=find(xX.iiS==i & xX.iiN<2); 
    X=xX.X(Sess(i).row,j); 
    y=Y(Sess(i).row,:);   % Get data from this block 
    beta(xX.iiS==i & xX.iiN==2,:)=mean(y,1); 
    ym=bsxfun(@minus,y,mean(y,1)); 
    for n=1:size(Y,2) 
        % Estimate hyperparameters 
        [V,h(:,n,i),k]=spmj_reml_sc2(ym(:,n)*ym(:,n)',{xX.Q{1,:}},1,8); 
                
        h(:,n,i)./xX.sh(i,:)'; 
        % Estimate beta weights 
        g=zeros(size(X,2),1);
        for hh=1:size(h,1)-1    
            g(xX.Gc{i,hh},1)=h(hh); 
        end;
        beta(j,n)=diag(g)*X'*(V\ym(:,n)); 
    end; 
end; 

% Check crossvalidation 
for i=1:nSess

    %  RSS(1,:)=RSS(1,:)+sum(ym.^2); 
    % RSS(2,:)=RSS(
    % 
    % keyboard;
    
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
varargout={beta,RSS,betaStats,mean(h,3)}; 