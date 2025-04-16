% Create the same test matrix in MATLAB
test_matrix = [
  4, 2, 1;
  2, 5, 3;
  1, 3, 6
];

% Test the get_sigma_lbm and demean functions
dist_mat = [
  0, 1, 2;
  1, 0, 3;
  2, 3, 0
];

function [sigma_dm] = demean(sigma)

  % With Var(X) = sigma .. comput the variance of Var(X-Xbar);
  n = size(sigma,1);
  sigma_dm = sigma-repmat(mean(sigma,1),n,1);
  sigma_dm = sigma_dm - repmat(mean(sigma_dm,2),1,n);

end

function [sigma_lbm] = get_sigma_lbm(distmat)
  % Compute the LBM covariance matrix from distmat using the first location as the origin
  n = size(distmat,1);
  sigma_lbm = 0.5*(repmat(distmat(:,1),1,n) + repmat(distmat(1,:),n,1) - distmat);
end


% Run your functions
sigma_lbm = get_sigma_lbm(dist_mat);
sigma_lbm_dm = demean(sigma_lbm);

% Then test Get_R
[R, DS] = Get_R(sigma_lbm_dm, 2);
disp('MATLAB values:');
disp(DS);
disp(R);

