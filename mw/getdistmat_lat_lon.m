function [distmat] = getdistmat_lat_lon(s)
% S = [lat lon]; get distance between points .. lat/lon points on sphere
% 
n = size(s,1);
distmat = NaN(n,n);
for j = 1:n
        tmp = repmat(s(j,:),n,1);
        distmat(:,j) = distance(s,tmp);
end
    
end

