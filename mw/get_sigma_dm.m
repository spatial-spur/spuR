function [sigma_dm] = get_sigma_dm(distmat,c)
% Computes demeaned version of sigma(c)
n = size(distmat,1);
sigma = exp(-c*distmat);
% Demean
sigma_dm = sigma-repmat(mean(sigma,1),n,1);
sigma_dm = sigma_dm - repmat(mean(sigma_dm,2),1,n);

end