function [eff,meanITI,diagM]=spmj_testdesign_param(ITI_type,varargin);
% Example function of how to find an optimal design for different
% questions
% Steps involved:
% Make a D-structure of Blocks, onsets and parmetric modulations 
% Assemble design into an SPM-structure: spmj_assemble 
% Make the X-matrix with:                spmj_fmri_design
% Calculate efficiency and may display:  spmj_efficiency, spmj_desMtx

type='normal';
random=1;
% Equal spaced 
switch(ITI_type)
    case 'equal'
        for i=1:100
            meanITI(i)=mod(i,10)+2;
            ITIs{i}=meanITI(i);
        end;
    case 'exponential'
        for i=1:100
            meanITI(i)=mod(i,10)+2;
            ITIs{i}=exponential_iti(meanITI(i),1,150);
        end;
    case 'miniblock'
        for i=1:100
            block_length=varargin{1};
            inITI=varargin{2};
            meanITI(i)=mod(i,10)+2;
            rest=meanITI(i)*(block_length)-inITI*(block_length-1);
            ITIs{i}=[repmat(inITI,block_length-1,1);rest];
            random=0;
        end;
end;

for i=1:length(ITIs);
    D=make_exp(ITIs{i},random,'blocked'); % generate the trial sequence and TR 
    % put into SPM.Session structure
    SPM=spmj_assemble(D,'nscan',[200 200],...
           'event','type',1,'move','P',...
            'basisfunc','hrf'); 
    SPM=spmj_fmri_design(SPM,'one_reg','demean','nonverbose'); % make the X-matrix
    ITIs{i};
    W=ones(400,1);
    for b=1:length(SPM.nscan)
        W(SPM.Sess(b).row(SPM.Sess(b).U(1).ons+1))=0.4;
    end;        
    W=diag(W);
    [eff(i),diagM(:,i)]=spmj_efficiency(SPM,0,'W',W);
end;




function D1=make_exp(ITIs,random,type,varargin)
if nargin<2
    random=1;
end;
ITI=mean(ITIs);
numtrials=floor(200/ITI);
if length(ITIs)<numtrials
    ITIs=repmat(ITIs,ceil(numtrials./length(ITIs)),1);
end;
D1=[];
for b=1:2
    if (random==1)
        ITI=sample_wor(ITIs,numtrials);
    else 
        ITI=ITIs(1:numtrials,1);
    end;
    D.TR=1;
    for i=2:numtrials
        D.TR(i,1)=D.TR(i-1,1)+ITI(i-1);
    end;
    D.type=ones(numtrials,1);
    % Now make parameteric modulation 
    switch (type)
        case 'random'
            D.P=normrnd(0,1,numtrials,1);
        case 'blocked'
            D.P=mod(b,2)*ones(numtrials,1);
        case 'event' 
            % Blocked trick trials 
            D.P=zeros(numtrials,1);
            block_length=varargin{1};
            a=sample_wor(1:floor(numtrials/block_length),max(1,floor(numtrials/(block_length*10))));
            for n=0:block_length-1
                D.P(a*block_length-n)=1;
            end;
    end;
    D.BN=ones(numtrials,1)*b;
    i=find(D.TR<200);
    D=getrow(D,i);
    D1=addstruct(D1,D);
end;
