function [cbar,cmax,W] = getcbar_W_large_n(rhobar,rhobar_min,s,latlongflag,qmax,m,N)
% Find cbar and W when n is large
 n = size(s,1);
 d = size(s,2);

 s_mN = NaN(m,d,N);
 distmat_mN = NaN(m,m,N); 
 r = s;
 for i = 1:N;
   r = jumble_s(r,m);
   s_tmp = r(1:m,:);
   distmat_tmp = getdistmat(s_tmp,latlongflag);
   s_mN(:,:,i) = s_tmp;
   distmat_mN(:,:,i) = distmat_tmp;
  end;
 
  cbar = getcbar_large_n(rhobar,distmat_mN);
  cmax = getcbar_large_n(rhobar_min,distmat_mN);
  
  % Compute Eigenvectors
  Wall = NaN(n,qmax*N); % save 1st evec as 1:N, then second as N+1:2N, ...
  for i = 1:N;
    s_tmp = squeeze(s_mN(:,:,i));
    distmat_tmp = squeeze(distmat_mN(:,:,i));
    W_tmp = getW(distmat_tmp,cbar,qmax);
    W_tmp = W_tmp(:,2:end);  % eliminate constant 
    % Nystrom Extension
    dist = zeros(n,m);
    for j = 1:d;
       dist = dist + (repmat(s(:,j),1,m)-repmat(s_tmp(:,j)',n,1)).^2;
    end;
    distmat = sqrt(dist);
    v = exp(-cbar*distmat);
    W = v*W_tmp;
    W = W-mean(W);
    tmp = sqrt(sum(W.^2));
    ii = (i:N:qmax*N);
    Wall(:,ii) = W./tmp;
   end;
   % Approximate Eigenspace
   r = NaN(n,qmax);
   for q = 1:qmax;
     [~,r_q,~]=pca(Wall(:,1:q*N),'NumComponents',1);  % Note PCA demeans by default
     r_q = r_q/sqrt(sum(r_q.^2));                     % Unit Length
     Wall = Wall - r_q*(r_q'*Wall);                   % Orthogonalize Wall wrt r_q
     r(:,q) = r_q;
   end;
   W = [ones(n,1)/sqrt(n) r];

end

