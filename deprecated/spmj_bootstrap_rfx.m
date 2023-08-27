function spmj_bootstrap_rfx(contrasts,punc,num_subj,varargin);
num_iter=1000;
mask=[];
threshold=tinv(1-punc,num_subj-1); % Uncorrected therhold to test

vararginoptions(varargin,{'num_iter','mask'});
T=[];

if (isempty(contrasts))
    contrasts=spm_select([1 inf],'image','select contrasts'); 
end; 
V=spm_vol(char(contrasts));
DATA=spm_read_vols(V);
G=size(DATA,4);

if (isempty(mask))
    mask=spm_select([1],'image','select mask'); 
end; 
    MV=spm_vol(mask);
    M=spm_read_vols(MV);


for n=1:num_iter
    a=sample_wr([1:G],num_subj);
    for i=1:num_subj
        X(:,:,:,i)=DATA(:,:,:,a(i)).*(unidrnd(2)*2-3).*M;
    end;
    m=mean(X,4);
    sd=std(X,0,4);
    t=m./(sd/sqrt(num_subj-1));
    indx=find(t>threshold);
    [x,y,z] = ind2sub(size(M),indx);
    D.Tmax=max(t(:));
    if isempty(x)
        D.numCl=0;
        D.maxSize=0;
    else
        A=spm_clusters([x y z]');
        vox=pivottable(A',[],A','length');
        
        D.numCl=length(vox);
        D.maxSize=max(vox);
    end;
    T=addstruct(T,D);
end;

fprintf('\nUncorrect threshold t(%d)=%2.3f  p=%2.3f\n',num_subj-1,threshold,punc);
fprintf('Corrected height-threshold (p<0.05):%2.3f\n',prctile(T.Tmax,95));
fprintf('Corrected Cluster-size threshold (at uncorrected height threshold):%2.3f\n',prctile(T.maxSize,95));


