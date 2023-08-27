function [S,R]=spmj_test_shear(M);
% Test for shear in transformation matrix 
R         = M(1:3,1:3);
vx        = sqrt(sum(M(1:3,1:3).^2));
vx(vx==0) = 1;
R         = R * diag(1./vx);

% Ensure that R is O(3)
[U,S,V] = svd(R);
R       = U*V';
if any(abs(diag(S)-1)>1e-3), warning('QFORM0 representation has been rounded.'); end;