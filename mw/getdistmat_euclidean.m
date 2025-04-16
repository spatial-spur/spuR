function [distmat] = getdistmat_euclidean(s)
% Get matrix of distances from location matrix
% 
n = size(s,1);
d = size(s,2);
dist = zeros(n,n);
for j = 1:d;
 dist = dist + (repmat(s(:,j),1,n)-repmat(s(:,j)',n,1)).^2;
end;
distmat = sqrt(dist);
    
end

