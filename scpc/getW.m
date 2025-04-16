function [W] = getW(distmat,c,qmax)
% Computes first qmax eigenvectors of demeaned Sigma(c0) and stores them in W, along with column vector of ones
% all columns of W are normalized to length 1
n = size(distmat,1);
sigma = exp(-c*distmat);
% Demean
sdm = sigma-repmat(mean(sigma,1),n,1);
sdm = sdm - repmat(mean(sdm,2),1,n);
% Eval
q = min([n;qmax]);
[V,D] = eig(sdm,'vector');
[D,ii] = sort(D,'descend');
DS = D(1:q);
r = V(:,ii(1:q));
W = [ones(n,1)/sqrt(n) r];

end