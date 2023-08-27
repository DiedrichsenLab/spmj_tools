function spmj_testdesign
% Example function of how to find an optimal design for different
% questions
% Steps involved:
% Make a D-structure of Blocks, onsets and parmetric modulations 
% Assemble design into an SPM-structure: spmj_assemble 
% Make the X-matrix with:                spmj_fmri_design
% Calculate efficiency and may display:  spmj_efficiency, spmj_desMtx

type='normal';
for i=[2:20];
    ITIs=i;
    D=make_exp(ITIs,type); % generate the trial sequence and TR 
    % put into SPM.Session structure
    SPM=spmj_assemble(D,'nscan',[200],...
           'event','type',1,'go','P',...
           'event','type',2,'nogo',{}); 
    SPM=spmj_fmri_design(SPM,'one_reg','demean'); % make the X-matrix
    [eff,diagM]=spmj_efficiency(SPM,1);
end;
keyboard;



function D=make_exp(ITIs,type)
ITI=mean(ITIs);
numtrials=floor(140/ITI);
numnogo=floor(numtrials*0);

D.TR=[1:5:200]';
D.BN=ones(length(D.TR),1);
D.type=unidrnd(3,40,1);
D.P=exp(1-D.TR/50);


