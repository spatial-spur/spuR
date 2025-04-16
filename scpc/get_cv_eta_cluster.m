function [cv_vec_eta] = get_cv_eta_cluster(x,Z,Xi_x,Xi_xs,xnorm_vec,r,cbar,cmax,distmat,alpha_vec)
 V = x;
 if size(Z,1) > 1
  V = [x Z];
 end
 n = size(r,1);
 q = size(r,2);
 W_eta_1 = xnorm_vec/sqrt(n);
 tmp = Xi_x*r;
 btmp = V\tmp;
 tmp = tmp-V*btmp;
 W_eta_2 = Xi_xs'*tmp;
 W_eta = [W_eta_1 W_eta_2];
 omega_mat_eta = form_omega_mat(cbar,cmax,q,W_eta,distmat);
 cv_vec_eta = NaN(length(alpha_vec),1);
 for ic = 1:length(alpha_vec)
   cv_vec_eta(ic)=findcv(omega_mat_eta,alpha_vec(ic),q);
 end
end