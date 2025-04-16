function plotrslt = c_ci_plot(Y,distmat,emat,rho_grid,rho_grid_ha)

% Construct confidence sets for c

q = size(emat,1);
n = size(distmat,1);
n_rho = length(rho_grid);
n_rho_ha = length(rho_grid_ha);

% BM covariance matrix (approximation for demeanded value)
rho_bm = 0.999;
c_bm = getcbar(rho_bm,distmat);
sigdm_bm = get_sigma_dm(distmat,c_bm);

% Construct R .. eigenvectors for low-frequency weights
[R lam]= Get_R(sigdm_bm,q);

% Step 1 .. compute Omega matrices for each value of rho_vec
ch_om_mat = NaN(q,q,n_rho);
ch_omi_mat = NaN(q,q,n_rho);
const_den_mat = NaN(n_rho,1);

for i = 1:n_rho
    rho = rho_grid(i);
    om = eye(q);
    if rho > 0
      c = getcbar(rho,distmat);
      sigdm = get_sigma_dm(distmat,c);
      om = R'*sigdm*R;
    end
    ch_om_mat(:,:,i) = chol(om);
    omi = inv(om);
    ch_omi_mat(:,:,i) = chol(omi);
    const_den_mat(i)=sqrt(det(omi))*0.5*gamma(q/2)/(pi^(q/2));
end

ch_omi_mat_ha = NaN(q,q,n_rho);
const_den_mat_ha = NaN(n_rho,1);
for i = 1:n_rho_ha
    rho = rho_grid_ha(i);
    om = eye(q);
    if rho > 0
      c = getcbar(rho,distmat);
      sigdm = get_sigma_dm(distmat,c);
      om = R'*sigdm*R;
    end
    omi = inv(om);
    ch_omi_mat_ha(:,:,i) = chol(omi);
    const_den_mat_ha(i)=sqrt(det(omi))*0.5*gamma(q/2)/(pi^(q/2));
end

% Construct critical values for each value of c_grid
sz_vec = [0.10 0.05 0.01]';
cv_mat = NaN(length(sz_vec),n_rho);
pv_mat = NaN(n_rho,1);
X = R'*Y;
den_ho_X_mat = NaN(n_rho,1);
den_ha_avg_X_mat = NaN(n_rho,1);
for i = 1:n_rho
    ch_null = squeeze(ch_om_mat(:,:,i));
    e = ch_null'*emat;
    ch_omi = squeeze(ch_omi_mat(:,:,i));
    const_den = const_den_mat(i);
    xc = ch_omi*e;
    den_ho = const_den*((sum(xc.^2,1)).^(-q/2))';
    den_ha_mat = NaN(size(emat,2),n_rho_ha);
    % Carry out for Data
    Xc = ch_omi*X;
    den_ho_X = const_den*((sum(Xc.^2,1)).^(-q/2))';
    den_ha_mat_X = NaN(n_rho_ha,1);
    for j = 1:n_rho_ha
        ch_omi = squeeze(ch_omi_mat_ha(:,:,j));
        const_den = const_den_mat_ha(j);
        xc = ch_omi*e;
        den_ha_mat(:,j) = const_den*((sum(xc.^2,1)).^(-q/2))';
        Xc = ch_omi*X;
        den_ha_mat_X(j,1) = const_den*((sum(Xc.^2,1)).^(-q/2));
    end
    den_ha_avg = mean(den_ha_mat,2);
    lr = den_ha_avg./den_ho;
    cvalue = prctile(lr,100*(1-sz_vec));
    cv_mat(:,i) = cvalue;
    den_ha_avg_X = mean(den_ha_mat_X);
    lr_X = den_ha_avg_X/den_ho_X;
    pv_mat(i)=mean(lr>lr_X);
    den_ha_avg_X_mat(i)=den_ha_avg_X;
    den_ho_X_mat(i) = den_ho_X;
end


% Save Variables for Plotting
ln_den = log(den_ho_X_mat);
cv_lines = repmat(den_ha_avg_X_mat,1,3)./cv_mat';
ln_cv_lines = log(cv_lines);
plotrslt = [ln_den pv_mat ln_cv_lines];

end