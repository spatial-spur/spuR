function [distmat] = getdistmat(s,latlongflag)
% Get matrix of distances from location matrix
% 

if latlongflag == 0
 n = size(s,1);
 d = size(s,2);
 dist = zeros(n,n);
 for j = 1:d;
  dist = dist + (repmat(s(:,j),1,n)-repmat(s(:,j)',n,1)).^2;
 end;
 distmat = sqrt(dist);
else
 n = size(s,1);
 distmat = NaN(n,n);
 for j = 1:n;
        tmp = repmat(s(j,:),n,1);
        distmat(:,j) = distance(s,tmp);
 end;
end

end

