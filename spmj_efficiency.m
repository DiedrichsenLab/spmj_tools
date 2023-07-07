function [eff,diagA,A]=spmj_efficiency(SPM,varargin)
% [eff,diagA,A]=spmj_efficiency(SPM,fig,varargin); 
% varargin: 
%   'W',W: weighting matrix (ACF structure)
%   'highpass',cuttoff: cutoff frequency of high-pass filter 
%   'C',C: contrast matrix 
%   'fig': Figure 
% Output: 
%   eff: Efficiency 1./trace(A)
%   diagA: Diagonal of A 
%   A: Variance-covariance matrix of contrasts (C'*var(beta)*C)
[nscan,np]=size(SPM.xX.X);
highpass=128;
SPM.xX.W= eye(nscan);
numreg=length(SPM.xX.iC);
numblock=length(SPM.xX.iB);
C=eye(numreg);
fig=0;

vararginoptions(varargin,{'W','highpass','C'},{'fig'});

% Make the high-pass-filters 
for b=1:length(SPM.Sess);
    K(b) = struct(	'HParam',	highpass,...
        'row',		SPM.Sess(b).row,...
        'RT',		SPM.xY.RT);
end;
SPM.xX.K = spm_filter(K);
SPM.xX.xKXs = spm_sp('Set',spm_filter(SPM.xX.K,SPM.xX.W*SPM.xX.X));	

% Add the zeros for the blocks and calculate the efficiency 
C=[C;zeros(numblock,size(C,2))];
if (fig>0)
    % O=spmj_desMtx(SPM.xX,'orthogonal','display');
    %keyboard;
    [eff,A]=spmj_desMtx(SPM.xX,'efficiency',C,'display');
    keyboard;
else 
    [eff,A]=spmj_desMtx(SPM.xX,'efficiency',C);
end;    
diagA=diag(A);
