function [rslt] = scpc(y_data,s,latlongflag,rhobar,ci_level)

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

n = size(y_data,1);
k = size(y_data,2);

q_qmax_tst = 0;

while q_qmax_tst == 0;
 qmax = min([n-1;qmax]);

 if n <= 2000

    % Compute matrix of distances
    distmat = getdistmat(s,latlongflag);

    % Determine cbar from locations a rhobar
    cbar = getcbar(rhobar,distmat);
    cmax = getcbar(rhobar_min,distmat);
  
    % Compute eigenvectors
    W = getW(distmat,cbar,qmax);
  
    % Get Grid of values of Omega
    omega_mat = form_omega_mat(cbar,cmax,qmax,W,distmat);
 
 else
   % Normalize locations
   [s,y_data] = normalize_s(s,y_data,latlongflag);
   
   % Compute cbar and W when n is large
   N = 20;
   m = 1000;
   M = 1000000; 
   [cbar,cmax,W] = getcbar_W_large_n(rhobar,rhobar_min,s,latlongflag,qmax,m,N);
   
   % Get Grid of values of Omega
   omega_mat = form_omega_mat_large_n(cbar,cmax,s,W,M);
 end

 % Find q and critical value
 [q,cv]=findqcv(omega_mat,1-ci_level);
 q
 cv
 
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
cv_vec = NaN(size(ci_level_vec,1),1);
for ic = 1:size(ci_level_vec,1);
  cv_vec(ic)=findcv(omega_mat,1-ci_level_vec(ic),q);
end;
    
    
 % Find Weights for PC
 r = W(:,2:q+1); % Weights for PC

 % Form SCPC_SE
 beta_hat = NaN(k,1);
 se_beta_hat = NaN(k,1);
 pv = NaN(k,1);
 ci = NaN(k,2);
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
   se_beta_hat(i) = sn;
   ts(i) = tstat;
   pv(i) = pvalue;
   ci(i,1) = ybar-cv*sn;
   ci(i,2) = ybar+cv*sn;
 end;

 rslt.beta_hat = beta_hat;
 rslt.se_beta_hat = se_beta_hat;
 rslt.tstat = ts;
 rslt.pvalue = pv;
 rslt.ci = ci;
 rslt.q = q;
 rslt.cv = cv;
 rslt.cv_vec = cv_vec;
 rslt.ci_level_vec = ci_level_vec;

end

