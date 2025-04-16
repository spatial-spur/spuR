function [sigma_lbm_dm] = get_sigma_lbm_dm(distmat)

% Compute the sigma_lbm with origin given by first location
sigma_lbm = get_sigma_lbm(distmat);

% Demean
sigma_lbm_dm = demean_sigma(sigma_lbm);

end