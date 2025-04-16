function [q,cv]=findqcv(omega_mat,alpha)
% Find critical value and optimal q
%   
global GQxw

qmax = size(omega_mat,1)-1;

% Step 1: Find critical value for each q = 1, 2, ... , qmax

% Find Critical Value for each q: 
% Note this is cv/sqrt(q) from paper. sqrt(q) is included at the bottom
% 
cvalue_vec = NaN(qmax,1);
for q = 1: qmax;
  cvalue_vec(q) = findcv(omega_mat,alpha,q); 
end;

% Step 3: Find expected length and minimize 
qvec = (1:1:qmax)';
Elength = (1./sqrt(qvec)).*cvalue_vec.*gamma((qvec+1)/2)./gamma(qvec/2);
[~,min_ind] = min(Elength);
q = qvec(min_ind);
cv = cvalue_vec(min_ind);

end

