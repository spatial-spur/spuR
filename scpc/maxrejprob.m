function [p,ii] = maxrejprob(cv,q,omega_mat)
% Find max rejection probability over a set of matrices in omega_mat
global GQxw
nc = size(omega_mat,3);
pvec = NaN(nc,1);
for i = 1:nc;
  A = squeeze(omega_mat(1:q+1,1:q+1,i));
  pvec(i) = rejprob(cv,A);
end;
[p,ii] = max(pvec);

end

