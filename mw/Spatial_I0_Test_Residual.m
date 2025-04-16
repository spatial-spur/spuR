function SP = Spatial_I0_Test_Residual(Y,X,distmat,emat)
q = size(emat,1);
n = size(distmat,1);

% Form M
M = eye(n)-X*(inv(X'*X))*X';

% BM covariance matrix (approximation for demeanded value)
rho_bm = 0.999;
c_bm = getcbar(rho_bm,distmat);
sigdm_bm = get_sigma_residual(distmat,c_bm,M);

% Construct R .. eigenvectors for low-frequency weights
[R lam]= Get_R(sigdm_bm,q);

% Value of om_ho for rhobar = rhobar_max
rho = 0.001;
c = getcbar(rho,distmat);
sigdm_rho = get_sigma_residual(distmat,c,M);
sigma = eye(n);
sigma_dm = sigma-repmat(mean(sigma,1),n,1);
sigma_dm = sigma_dm - repmat(mean(sigma_dm,2),1,n);
sigdm_wn = sigma_dm;

% Matrices used in analysis;
om_rho = R'*sigdm_rho*R;
om_wn = R'*sigdm_wn*R;
om_bm = R'*sigdm_bm*R;

% Find value of Ha parameter that yields (roughly 50% power)
om_i0 = om_rho;
om_ho = om_rho;
ha_parm = get_ha_parm_I0(om_ho,om_i0,om_bm,emat); % Note .. no need to change this for residual-based test
om_ha = om_i0+ha_parm*om_bm;

% Matrices for carrying out tests
ch_omi_ho = chol(inv(om_ho));
ch_omi_ha = chol(inv(om_ha));

% Get LR for data
n_y = size(Y,2);
LR = NaN(n_y,1);
for i = 1:n_y
  X = Y(:,i) - mean(Y(:,i));
  P = R'*X;
  y_P_ho = ch_omi_ho*P;
  y_P_ha = ch_omi_ha*P;
  q_P_ho = sum(y_P_ho.^2,1)';
  q_P_ha = sum(y_P_ha.^2,1)';
  LR(i) = q_P_ho./q_P_ha;
end

% Construct om_ho for a grid of values of rho
rho_min = 0.0001;
rho_max = 0.03;
n_rho = 30;
rho_grid = linspace(rho_min,rho_max,n_rho)';
ch_om_ho_mat = NaN(q,q,n_rho);
for i = 1:n_rho
    om_ho = eye(q);
    rho = rho_grid(i);
    if rho > 0
      c = getcbar(rho,distmat);
      sigdm_ho = get_sigma_residual(distmat,c,M);
      om_ho = R'*sigdm_ho*R;
    end
    ch_om_ho_mat(:,:,i)=chol(om_ho);
end

pvalue_mat = NaN(n_rho,n_y);
sz_vec = [0.01 0.05 0.10];
cvalue_mat = NaN(n_rho,3);

for ir = 1:n_rho
    ch_om_ho = squeeze(ch_om_ho_mat(:,:,ir));
    y_ho = ch_om_ho'*emat;
    y_ho_ho = ch_omi_ho*y_ho;
    y_ho_ha = ch_omi_ha*y_ho;
    q_ho_ho = sum(y_ho_ho.^2,1)';
    q_ho_ha = sum(y_ho_ha.^2,1)';
    lr_ho = q_ho_ho./q_ho_ha;
    cv_vec = prctile(lr_ho,100*(1-sz_vec));
    cvalue_mat(ir,:)=cv_vec;
    for j = 1:n_y
        pvalue_mat(ir,j)=mean(lr_ho>LR(j));
    end
end

cvalue = max(cvalue_mat)';
pvalue = max(pvalue_mat)';


SP.LR = LR;
SP.pvalue = pvalue;
SP.cvalue = cvalue;
SP.ha_parm = ha_parm;
SP.cvalue_mat = cvalue_mat;
SP.pvalue_mat = pvalue_mat;
SP.rho_grid = rho_grid;

end