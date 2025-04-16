function [omega_mat] = form_omega_mat_large_n(cbar,cmax,s,W,M)

 n = size(s,1);
 qmax = size(W,2)-1;
 % Compute Distance Matrix for a random subset of 
 % Draw indices for locations
 inds = raninds(n,M);
 dsM = s(inds(1:M),:)-s(inds(2:M+1),:);
 distvec = sqrt(sum(dsM.^2,2));  % Distance vector .. pairwise j and j-1
   
 % Step 2: Compute Omega Matrices for grid of values of c >= cbar
 nc = ceil(log(cmax/cbar)/log(1.2));
 nc = max([2;nc]);
 c_grid = cbar*(1.2.^(0:1:nc)');
 omega_mat = NaN(qmax+1,qmax+1,nc);
 for i = 1:nc;
      tmp1 = exp(-c_grid(i)*distvec); 
      tmp2 = repmat(tmp1,1,qmax+1).*W(inds(2:M+1),:);
      tmp = W(inds(1:M),:)'*tmp2;
      omega_mat(:,:,i) = eye(qmax+1) + (n*(n-1)/(2*M))*(tmp+tmp');
 end;

end

