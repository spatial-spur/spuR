function [rslt] = cscpc(Y,X,s,latlongflag,rhobar,ci_level)

n_big = 3500;

n = size(X,1);
k = size(X,2);

if n > n_big
  fprintf(['n = ' num2str(n) ' is larger than threshold n = ' num2str(n_big) ':  Using large n method \n'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (a) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n > n_big
 % Normalize locations .. useful for invariance for some probabilistic
 % calculations with random selection of locations
 [s,tmp] = normalize_s(s,[Y X],latlongflag);
 Y = tmp(:,1);
 X = tmp(:,2:end);
end

%%%%%%%%%%% (b) Estimate Regression and compute HR standard errors %%%%%%%%%%%%%%
% Form XX, XY, etc;
XX = (X'*X);
XY = (X'*Y);
XX_inv = inv(XX);
S_XX = (X'*X)/n;
S_XX_inv = inv(S_XX);

% OLS
beta_hat = XX_inv*XY;

% Residuals
u_hat = Y-X*beta_hat;
Xu = X.*repmat(u_hat,1,k);  % X*u

% White Covariance
Xu_Xu = Xu'*Xu;
V_beta = XX_inv*Xu_Xu*XX_inv';
% STATA DF adjustment
V_beta = V_beta*n/(n-k);
se_beta_hat = sqrt(diag(V_beta));

% Save Results
rslt.beta_hat = beta_hat;           % Beta_hat
rslt.se_hr_beta_hat = se_beta_hat;  % HR Beta_hat

%%%%%%%%%%%%%%%%%%%%% (c) SCPC, R, q, and critical values %%%%%%%%%%%%%%%%%%
y_data = repmat(beta_hat',n,1)+Xu*S_XX_inv';     % j'th column is beta_hat_j+appropiate linear combination of Xu

nobs_s = size(s,1);
s1 = unique(s,'rows');
nobs_unique = size(s1,1);

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

n = size(y_data,1);
k = size(y_data,2);

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
   cv_vec_eta = get_cv_eta(x,r,Z,cbar,cmax,distmat,1-[ci_level;ci_level_vec]);
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
 
 for i = 1:k;
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
 end;

 rslt.beta_hat = beta_hat;
 rslt.se_beta_hat_scpc = se_beta_hat_scpc;
 rslt.tstat = ts;
 rslt.q = q;
 rslt.ci_level_vec = ci_level_vec;
 rslt.pvalue_scpc = pv_scpc;
 rslt.ci_scpc = ci_scpc;
 rslt.cv_scpc = cv_scpc;
 rslt.cv_vec_scpc = cv_vec_scpc;
 rslt.ci_cscpc = ci_cscpc;
 rslt.cv_cscpc = cv_cscpc;
 rslt.cv_vec_cscpc = cv_vec_cscpc;
 rslt.n = n;
 rslt.n_big = n_big;


end

