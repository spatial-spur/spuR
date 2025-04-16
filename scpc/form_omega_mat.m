function [omega_mat] = form_omega_mat(cbar,cmax,qmax,W,distmat)

% Step 1: Compute Omega Matrices for grid of values of c >= cbar
nc = ceil(log(cmax/cbar)/log(1.2));
nc = max([2;nc]);
c_grid = cbar*(1.2.^(0:1:nc)');
omega_mat = NaN(qmax+1,qmax+1,nc);
for i = 1:nc;
    omega_mat(:,:,i) = W'*(exp(-c_grid(i)*distmat))*W;
end;

end
