function [cv_vec_eta] = get_cv_eta(x,r,Z,cbar,cmax,distmat,alpha_vec)
 V = x;
 if size(Z,1) > 1
  V = [x Z];
 end
 n = size(r,1);
 q = size(r,2);
 sgn_x = (x >= 0) - (x < 0);
 W_eta_1 = abs(x)/sqrt(n);
 tmp = (repmat(x,1,q).*r);
 btmp = V\tmp;
 tmp = tmp-V*btmp;
 W_eta_2 = repmat(sgn_x,1,q).*tmp;
 W_eta = [W_eta_1 W_eta_2]; 
 omega_mat_eta = form_omega_mat(cbar,cmax,q,W_eta,distmat);
 cv_vec_eta = NaN(length(alpha_vec),1);
 for ic = 1:length(alpha_vec)
   cv_vec_eta(ic)=findcv(omega_mat_eta,alpha_vec(ic),q);
 end
end

