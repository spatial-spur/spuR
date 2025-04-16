function [X] = sqrt_psd(U)
% U is a PSD Matrix.  U = X*X';  
% Computed using eigenvalues (probably not efficient)
%
% Use Cholesky if possible, if not use eigenvector
[Y,flag]=chol(U);
if flag == 0
  X = Y';
else
 [V,D]=eig(U);
 d = diag(D);
 ii = d<0;
 d(ii==1) = 0;
 D1 = diag(sqrt(d));
 V = real(V);
 X = V*D1;
end

end

