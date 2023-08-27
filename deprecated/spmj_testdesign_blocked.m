function E=spmj_testdesign_blocked(varargin);
% Example function of how to find an optimal design for different
% questions
% Steps involved:
% Make a D-structure of Blocks, onsets and parmetric modulations 
% Assemble design into an SPM-structure: spmj_assemble 
% Make the X-matrix with:                spmj_fmri_design
% Calculate efficiency and may display:  spmj_efficiency, spmj_desMtx

block_length=[2:20];    % 10 to 40 s
gap_length=[2:20];      % 10 to 40 s
num_cond=4;
N=5;                   % number of simulations per condition
nscan=[160 160];
fig=0;
C=[eye(4);1 -1 0 0;1 0 -1 0;1 0 0 -1;0 1 -1 0;0 1 0 -1;0 0 1 -1]';
if (nargin>0)
    PARAM=struct(varargin{:});
    if exist('PARAM', 'var') & ~isempty(PARAM)
        fnames = fieldnames(PARAM);
        for fnum = 1:length(fnames);
            eval([fnames{fnum} ' = getfield(PARAM, fnames{fnum});']);
        end;
    end;
end;

E.eff=[];
E.diagM=[];E.block_length=[];E.gap_length=[];
for i=1:length(block_length);
    for n=1:N
        D=make_exp(block_length(i),gap_length(i),num_cond,nscan);
        SPM=spmj_assemble(D,'nscan',nscan,...
            'event','con1','type',1,[],...
            'event','con2','type',2,[],...
            'event','con3','type',3,[],...
            'event','con4','type',4,[],...
            'basisfunc','hrf'); 
        SPM=spmj_fmri_design(SPM,'one_reg','demean','nonverbose'); % make the X-matrix
        [E.eff(end+1,1),E.diagM(:,end+1)]=spmj_efficiency(SPM,fig,'C',C,'high_pass',200);
        E.block_length(end+1,1)=block_length(i);
        E.gap_length(end+1,1)=gap_length(i);
    end;
end;

function D1=make_exp(block_length,gap_length,num_cond,nscan)
D1=[];
for b=1:length(nscan)
    num_blocks=floor(nscan(b)/(block_length+gap_length));
    D.TR=2;
    for i=2:num_blocks
        D.TR(i,1)=D.TR(i-1,1)+(block_length+gap_length);
    end;
    samp=repmat([1:num_cond]',ceil(num_blocks/num_cond),1);
    samp=samp(1:num_blocks);
    D.type=sample_wor(samp,num_blocks,1);
    D.BN=ones(num_blocks,1)*b;
    D.dur=ones(num_blocks,1)*block_length;
    D1=addstruct(D1,D);
end;
