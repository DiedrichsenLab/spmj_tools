function varargout=spmj_checkFilter(what,varargin); 
% Investigates the influence of temporal filtering on first-level GLM 

switch (what)
    case 'getW' 
        SPM=varargin{1}; 
        nSess= length(SPM.Sess) 
        for n=1:nSess 
            indxR = SPM.Sess(n).row; 
            indxC = SPM.Sess(n).col; 
            V(:,:,n)=full(SPM.xVi.V(indxR,indxR)); 
            Vx(:,:,n)=full(SPM.xX.V(indxR,indxR)); 
            W(:,:,n)=real(((V(:,:,n)+V(:,:,n)')/2)^(-1/2)); 
            Wx(:,:,n)=full(SPM.xX.W(indxR,indxR)); 
        end;
        if nargout==0
        subplot(2,2,1); 
        plot(squeeze(V(350,:,:))); 
        subplot(2,2,2); 
        plot(squeeze(Vx(350,:,:))); 
        subplot(2,2,3); 
        plot(squeeze(W(350,:,:))); 
        subplot(2,2,4); 
        plot(squeeze(Wx(350,:,:))); 
        end; 
        varargout={Wx}; 
    case 'checkX' 
        SPM=varargin{1}; 
        indxR = SPM.Sess(1).row; 
        indxC = SPM.Sess(1).col; 
        XX = full(SPM.xX.X(indxR,indxC)); 
        % XX = bsxfun(@minus,XX,mean(XX)); 
        X = [XX ones(length(indxR),1)]; 
        W = spmj_checkFilter('getW',SPM); 
        WX=W(:,:,1)*X; 
        plot(sum(WX(:,1:end-1),2)); 
        hold on; 
        plot(sum(WX(:,end),2));
        hold off; 
        
        B=ones(21,1);
        Xp = X; 
        Xp(X<0)=0; 
        Y=Xp*B; 
        
        B_hat1 = pinv(X)*Y;  
        B_hat2 = pinv(WX)*W(:,:,1)*Y; 
        
        keyboard; 
    case 'checkReg' 
        SPM = varargin{1}; 
        Y   = varargin{2}; 
        W = spmj_checkFilter('getW',SPM); 
        for i=1:8 
            indxR = SPM.Sess(i).row; 
            indxC = SPM.Sess(i).col; 
            XX = full(SPM.xX.X(indxR,indxC)); 
            X = [XX ones(length(indxR),1)]; 
            B1(:,:,i) = pinv(W(:,:,i)*X)*W(:,:,i)*Y(indxR,:); 
            B2(:,:,i) = pinv(X)*Y(indxR,:); 
        end; 
        keyboard; 
end; 