function [eff,meanITI,diagM]=spmj_testdesign_dartmouth(varargin);
% Example function of how to find an optimal design for different
% questions
% Steps involved:
% Make a D-structure of Blocks, onsets and parmetric modulations 
% Assemble design into an SPM-structure: spmj_assemble 
% Make the X-matrix with:                spmj_fmri_design
% Calculate efficiency and may display:  spmj_efficiency, spmj_desMtx
T=252; % Length of scans 

ITIs=[0.5:0.5:10]';
ITIs=repmat(ITIs,10,1);

% Contrast matrix: detect, block, unibi, congruent_incongurent
C=[1 1 1 1 1 1 1 1;1 1 1 1 -1 -1 -1 -1;1 1 -1 -1 1 1 -1 -1;0 0 1 -1 0 0 1 -1];

for i=1:length(ITIs)
    D=make_exp(ITIs(i),T); % generate the trial sequence and TR 
    % put into SPM.Session structure
    SPM=spmj_assemble(D,'nscan',floor(T/2.5),...
           'event','type',1,'uniLDir',[],...
           'event','type',2,'uniRDir',[],...
           'event','type',3,'biCoDir',[],...
           'event','type',4,'biInDir',[],...
           'event','type',5,'uniLSym',[],...
           'event','type',6,'uniRSym',[],...
           'event','type',7,'biCoSym',[],...
           'event','type',8,'biInSym',[],...
           'TR',2.5,...
            'basisfunc','hrf'); 
    SPM=spmj_fmri_design(SPM,'one_reg','demean','nonverbose'); % make the X-matrix
    [eff(i),diagM(:,i),A(:,:,i)]=spmj_efficiency(SPM,1);
    eff_con(:,i)=1./diag(C*A(:,:,i)*C');   
    keyboard; % take out if you run batch mode
end;
E.ITI=ITIs;
E.detect=eff_con(1,:)';
E.block=eff_con(2,:)';
E.unibi=eff_con(3,:)';
E.coinc=eff_con(4,:)';

keyboard;


function D=make_exp(ITI,T)
numtrials=floor(T/ITI);
numblocks=floor(numtrials/12)+1;
D.type=[];
TRIALS=[1 1 2 2 3 3 4 4 0 0 0 0];
for b=1:numblocks
    t=sample_wor(TRIALS,12);
    t1=t+mod(b,2)*4;t1(t==0)=0;
    D.type(end+1:end+12,1)=t1;
end;
D.sec=[0:length(D.type)-1]';
D.sec=D.sec*ITI;
i=find(D.sec<T);
D=getrow(D,i);