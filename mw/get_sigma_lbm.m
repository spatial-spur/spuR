function [sigma_lbm] = get_sigma_lbm(distmat)
    % Compute the LBM covariance matrix from distmat using the first location as the origin
    n = size(distmat,1);
    sigma_lbm = 0.5*(repmat(distmat(:,1),1,n) + repmat(distmat(1,:),n,1) - distmat);
end
