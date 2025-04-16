function SP = Spatial_I1_Test_Residual(Y,X,distmat,emat)
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

% Matrices used in analysis;
om_ho = R'*sigdm_bm*R;

% Find value of Ha parameter that yields (roughly 50% power)
ha_parm = get_ha_parm_I1_residual(om_ho,distmat,R,emat,M);
sigdm_ha = get_sigma_residual(distmat,ha_parm,M);
om_ha = R'*sigdm_ha*R;

% Get Draws of LR under ho
ch_om_ho = chol(om_ho);
omi_ho = inv(om_ho);
omi_ha = inv(om_ha);
ch_omi_ho = chol(omi_ho);
ch_omi_ha = chol(omi_ha);
y_ho = ch_om_ho'*emat;
y_ho_ho = ch_omi_ho*y_ho;
y_ho_ha = ch_omi_ha*y_ho;
q_ho_ho = sum(y_ho_ho.^2,1)';
q_ho_ha = sum(y_ho_ha.^2,1)';
lr_ho = q_ho_ho./q_ho_ha;

sz_vec = [0.01 0.05 0.10];
cv_vec = prctile(lr_ho,100*(1-sz_vec));

n_y = size(Y,2);
LR = NaN(n_y,1);
pvalue = NaN(n_y,1);
for i = 1:n_y
  X = Y(:,i) - mean(Y(:,i));
  P = R'*X;
  LR(i) = (P'*omi_ho*P)/(P'*omi_ha*P);
  pvalue(i) = mean(lr_ho>LR(i));
end

SP.LR = LR;
SP.pvalue = pvalue;
SP.cv_vec = cv_vec;
SP.ha_parm = ha_parm;

end