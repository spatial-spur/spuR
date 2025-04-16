function [distmat] = getdistmat_normalized(s,latlongflag)
% Get matrix of distances from location matrix
% Normalize the matrix so the largest distance is 1

if latlongflag == 0
   % Use Euclidean distance
   distmat = getdistmat_euclidean(s);
else
   % Use Latitute/Longitude distance .. distance on a sphere
   distmat = getdistmat_lat_lon(s);
end
 
% Normalize the matrix
distmat = distmat/max(max(distmat));      

end

