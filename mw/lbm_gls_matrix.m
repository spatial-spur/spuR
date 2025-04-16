function [LBMGLS_mat]=lbm_gls_matrix(s,latlongflag)
% Compute the lbm_gls matrix 
% Input:
%   s: (n,2) matrix of locations
%   latlongflag: 1 if the input is latitute and longitude, 0 if euclidean distance is to be used

small = 1.0e-10;

% Compute distance matrix 
distmat = getdistmat_normalized(s,latlongflag);
% BM covariance matrix demeanded
sigma_lbm_dm = get_sigma_lbm_dm(distmat);

% Compute the eigenvectors and eigenvalues 
[V,d] = eig(sigma_lbm_dm,'vector');
[eval, ii] = sort(d,'descend');
evec = V(:,ii);

% Compute the GLS matrix
ii = eval > small;
eval = eval(ii==1);
evec = evec(:,ii==1);
dsi = 1./sqrt(eval);
Dsi = diag(dsi);
LBMGLS_mat = evec*Dsi*evec';

end