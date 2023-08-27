function spmj_testdesign_evoked
% Example function of how to find an optimal design for different
% questions
% Steps involved:
% Make a D-structure of Blocks, onsets and parmetric modulations 
% Assemble design into an SPM-structure: spmj_assemble 
% Make the X-matrix with:                spmj_fmri_design
% Calculate efficiency and may display:  spmj_efficiency, spmj_desMtx

type='normal';
ITI1=exponential_iti(2,1,150);
ITI2=exponential_iti(3,1,150);
ITI3=exponential_iti(4,1,150);
ITI4=exponential_iti(5,1,150);
ITI5=exponential_iti(5,1,150);
ITI6=exponential_iti(5,1,150);
ITI7=exponential_iti(5,1,150);


ITIs={10,11,12,13,14,15,16,17,18,19,20,ITI2,ITI3,ITI4,ITI5,ITI6,ITI7};
for i=1:length(ITIs);
    D=make_exp(ITIs{i},type); % generate the trial sequence and TR 
    % put into SPM.Session structure
    SPM=spmj_assemble(D,'nscan',[200],...
           'event','type',1,'even',{},...
            'basisfunc','fourier'); 
    SPM=spmj_fmri_design(SPM,'one_reg','demean'); % make the X-matrix
    [eff(i),diagM(:,i)]=spmj_efficiency(SPM,1);
end;
keyboard;



function D=make_exp(ITIs,type)
ITI=mean(ITIs);
numtrials=floor(200/ITI);
if length(ITIs)<numtrials
    ITIs=repmat(ITIs,ceil(numtrials./length(ITIs)),1);
end;
ITI=sample_wor(ITIs,numtrials);
D.TR=1;
for i=2:numtrials
    D.TR(i,1)=D.TR(i-1,1)+ITI(i);
end;
D.type=ones(numtrials,1);
D.P=exp(1-D.TR/50);
D.BN=ones(numtrials,1);


