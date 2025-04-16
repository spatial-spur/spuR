function [sigma_dm] = get_sigma_residual(distmat,c,M)
% Computes demeaned version of sigma(c)
n = size(distmat,1);
sigma = exp(-c*distmat);
% Residual Version
sigma_dm = M*sigma*M';

end