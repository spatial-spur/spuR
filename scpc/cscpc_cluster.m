function [rslt] = cscpc_cluster(Y,X,Cid,S,latlongflag,feflag,rhobar,ci_level)

small = 1.0e-10;
n_big = 3500;

N = length(Y);
d = size(S,2);

% Get Locations .. and normalize
cid = unique(Cid,'rows');
n = length(cid);           % Number of clusters
s = NaN(n,d);
for i = 1:n;
  ii = Cid == cid(i);
  stmp = S(ii==1,:);
  s(i,:)= stmp(1,:);
end

if n > n_big
  fprintf(['n = ' num2str(n) ' is larger than threshold n = ' num2str(n_big) ':  Using large n method \n'])
end

if n > n_big
 % Normalize locations .. useful for invariance for some probabilistic
 % calculations with random selection of locations
 [s,cid] = normalize_s(s,cid,latlongflag);  % Note the rows of s and cid have been rearranged
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up variables for clusters, locations, etc. %%%%%%%%%%%%%%%%%%%%%%%
m_vec = NaN(n,1);          % Vector with i'th entry = m(i)
for i = 1:n;
    ii = Cid == cid(i);
    m = sum(ii);
    ytmp = Y(ii==1,:);
    xtmp = X(ii==1,:);
    ytmp_dm = ytmp-mean(ytmp);
    xtmp_dm = xtmp-mean(xtmp);
    if i == 1
        x = xtmp;
        y = ytmp;
        ydm = ytmp_dm;
        xdm = xtmp_dm;
    else
        x = [x;xtmp];
        y = [y;ytmp];
        xdm = [xdm;xtmp_dm];
        ydm = [ydm;ytmp_dm];
    end
    m_vec(i)= m;
end
if feflag == 0
    Y = y;
    X = [x ones(N,1)];  % Append constant as last regressor
elseif feflag == 1
    Y  = ydm;
    X = xdm;   % entity demeaned
end
k = size(X,2);

%%%%%%%%%%% (b) Estimate Regression and compute clustered standard errors %%%%%%%%%%%%%%
beta_hat = X\Y;    
e_hat = Y - X*beta_hat;
XX_inv = inv(X'*X);
S_XX_inv = n*XX_inv;
Xu = NaN(n,k);
for l = 1:n
    xl = get_xl_X(X,l,m_vec);
    el = get_xl_X(e_hat,l,m_vec);
    Xu(l,:)=el'*xl;
end

% Clustered Covariance
Suu = Xu'*Xu;
V_beta = XX_inv*Suu*XX_inv';
% STATA DF adjustment
V_beta = V_beta*(n/(n-1))*((N-1)/(N-k));
se_beta_hat_cluster = sqrt(diag(V_beta));

rslt.beta_hat = beta_hat;                        % Beta_hat
rslt.se_beta_hat_cluster = se_beta_hat_cluster;  % HR Beta_hat


%%%%%%%%%%%%%%%%%%%%% (c) SCPC, R, q, and critical values %%%%%%%%%%%%%%%%%%
y_data = repmat(beta_hat',n,1)+Xu*S_XX_inv';     % j'th column is beta_hat_j+appropiate linear combination of Xu = X'e by cluster

if rhobar < 0.001
   fprintf('rhobar = %6.4f \n',rhobar);
   error('rhobar cannot be less than 0.001');
elseif rhobar > 0.99
   fprintf('rhobar = %6.4f \n',rhobar);
   error('rhobar cannot be greater than 0.99');
end

% global variables
global GQxw;
GQxw = setGQxw;
global random_t;
random_t = 1;

rhobar_min = 0.00001;     % Smallest value of rho .. for calibrating largest c values

% Set qmax;
if rhobar >= 0.05
   qmax = 10;
 elseif rhobar >= 0.01
     qmax = 20;
 elseif rhobar >= 0.005
     qmax = 60
 else
     qmax = 120;
end


% Compute matrix of distances
distmat = getdistmat(s,latlongflag);

q_qmax_tst = 0;
while q_qmax_tst == 0
 qmax = min([n-1;qmax]);

 if n <= n_big

    % Determine cbar from locations a rhobar
    cbar = getcbar(rhobar,distmat);
    cmax = getcbar(rhobar_min,distmat);
    
    % Compute eigenvectors
    W = getW(distmat,cbar,qmax);
    

    % Get Grid of values of Omega
    omega_mat = form_omega_mat(cbar,cmax,qmax,W,distmat);

 else
   
   % Compute cbar and W when n is large
   N1 = 20;
   m1 = 1000;
   M1 = 1000000;
   [cbar,cmax,W] = getcbar_W_large_n(rhobar,rhobar_min,s,latlongflag,qmax,m1,N1);
   
   % Get Grid of values of Omega
   omega_mat = form_omega_mat_large_n(cbar,cmax,s,W,M1);
 
 end

 % Find q and critical value
 [q,cv_scpc]=findqcv(omega_mat,1-ci_level);
 

 if q == qmax
   if q == n-1
     q_qmax_tst = 1;
   else
     qmax = round(1.5*qmax);
   end
   qmax
 else
  q_qmax_tst = 1;
 end
 
end

% Find Critical values for other coverage rates .. hold q constant
ci_level_vec = [0.68;0.90;0.95;0.99];
cv_vec_scpc = NaN(size(ci_level_vec,1),1);
for ic = 1:size(ci_level_vec,1);
  cv_vec_scpc(ic)=findcv(omega_mat,1-ci_level_vec(ic),q);
end;
        
 % Find Weights for PC
 r = W(:,2:q+1); % Weights for PC


%%%%%%%%%%%%%%%%%%%%% (d) C-SCPC critical values %%%%%%%%%%%%%%%%%%

 % Construct CSCPC critical value for each regressor
 Xsave = X;  % Save X 
 Z = nan;
 cv_cscpc = NaN(k,1);
 cv_vec_cscpc = NaN(size(ci_level_vec,1),k);

 for ik = 1:k
   X = Xsave(:,ik);
   if k > 1
     Z = Xsave;
     Z(:,ik)=[];
   end
   x = X;
   if size(Z,1) > 1
     bZ = Z\X;
     x = X-Z*bZ;
   end
   xs = NaN(N,1);
   xnorm_vec = NaN(n,1);
   for l = 1:n
       xl = get_xl_X(x,l,m_vec);
       xnorm = sqrt(xl'*xl);
       xls = zeros(size(xl));
       if xnorm > small
          xls = xl/xnorm;
       end
       xs = put_xl_X(xs,xls,l,m_vec);
       xnorm_vec(l) = xnorm;
   end
   Xi_x = get_Xi(x,m_vec);
   Xi_xs = get_Xi(xs,m_vec);
   cv_vec_eta = get_cv_eta_cluster(x,Z,Xi_x,Xi_xs,xnorm_vec,r,cbar,cmax,distmat,1-[ci_level;ci_level_vec]);
   cv_cscpc(ik) = max([cv_scpc cv_vec_eta(1)]);
   for ic = 1:size(ci_level_vec,1)
      cv_eta = cv_vec_eta(ic+1);
      cv_vec_cscpc(ic,ik)=max([cv_vec_scpc(ic) cv_eta]);
   end 
 end 
   

 %%%%%%%%%%%%%%%%%%%% (e) t-stats and confidence intervals %%%%%%%%%%%%%%%%  
     
 % Form SCPC and CSCPC CI, etc.
 beta_hat = NaN(k,1);
 se_beta_hat_scpc = NaN(k,1);
 pv_scpc = NaN(k,1);
 ci_scpc = NaN(k,2);
 ci_cscpc = NaN(k,2);
 ts = NaN(k,1);
 
 for i = 1:k
   y = y_data(:,i);
   ybar = mean(y);
   uhat = y - ybar;
   pc = r'*uhat;
   shat = sqrt(sum(pc.^2)/q);
   sn = shat/sqrt(n);
   tstat = ybar/sn;
   pvalue = maxrejprob(abs(tstat/sqrt(q)),q,omega_mat);
   beta_hat(i) = ybar;
   se_beta_hat_scpc(i) = sn;
   ts(i) = tstat;
   pv_scpc(i) = pvalue;
   ci_scpc(i,1) = ybar-cv_scpc*sn;
   ci_scpc(i,2) = ybar+cv_scpc*sn;
   ci_cscpc(i,1) = ybar-cv_cscpc(i)*sn;
   ci_cscpc(i,2) = ybar+cv_cscpc(i)*sn;
 end

 rslt.beta_hat = beta_hat;
 rslt.se_beta_hat_scpc = se_beta_hat_scpc;
 rslt.tstat = ts;
 rslt.q = q;
 rslt.n = n;
 rslt.N = N;
 rslt.k = k;
 rslt.n_big = n_big;
 rslt.ci_level_vec = ci_level_vec;
 rslt.pvalue_scpc = pv_scpc;
 rslt.ci_scpc = ci_scpc;
 rslt.cv_scpc = cv_scpc;
 rslt.cv_vec_scpc = cv_vec_scpc;
 rslt.ci_cscpc = ci_cscpc;
 rslt.cv_cscpc = cv_cscpc;
 rslt.cv_vec_cscpc = cv_vec_cscpc;


end