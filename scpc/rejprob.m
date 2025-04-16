function [p] = rejprob(cv,A)
% Find max rejection probability over a set of matrices in omega_mat
global GQxw
  A(2:end,:) = -(cv^2)*A(2:end,:);
  w = eig(A);
  w = sort(w,'descend');
  eta = -w(2:end)/w(1);
  p = GQp(eta,GQxw);
end

