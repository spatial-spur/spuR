function [sigma_dm] = demean(sigma)

  % With Var(X) = sigma .. comput the variance of Var(X-Xbar);
  n = size(sigma,1);
  sigma_dm = sigma-repmat(mean(sigma,1),n,1);
  sigma_dm = sigma_dm - repmat(mean(sigma_dm,2),1,n);

end