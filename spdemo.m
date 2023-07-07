function sp = spdemo(arg1, arg2)
% Comparison of sparse vs full matrices
if (nargin ~= 2)
    error('Give order and density');
    return
end
S = sprandn(arg1,arg1,arg2);
F = full(S);
% Compare speed.

t0=cputime;
B = S * S;
stime=cputime-t0;

t0=cputime;
C = F * F;
ftime=cputime-t0;

sprintf('Order %d matrix, density %f: Full=%f, Sparse=%f diff=%f', arg1, ...
	arg2, ftime, stime,ftime-stime)