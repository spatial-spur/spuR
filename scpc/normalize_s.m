function [s,y] = normalize_s(s,y,latlongflag)

% normalizes locations to ensure invariance given probabilistic algorithm
% if locations are in lat/long format, just normalize longitudes (so invariance of rotations around polar axis)

if latlongflag == 0
  % normalize locations
   [~,s,~]=pca(s);  % Note PCA demeans by default
  % normalize sign
   if max(s(:,1)) ~= max(abs(s(:,1)))
       s = -s;
   end 
  % sort locations
   [~,ii]=sort(s(:,1));
   s = s(ii,:);
   y = y(ii,:);
   % put s in [0,1]^d
   s = s-min(s);
   s = s./max(max(s));  
end

if latlongflag == 1
   lat = s(:,1);
   lon = s(:,2);
   lon = lon-mean(lon);
   lon = mod(lon+180,360)-180;
   [~,ii] = sort(lon);
   lat = lat(ii);
   lon = lon(ii);
   y = y(ii,:);
   s = [lat lon];
end
   
  
end

