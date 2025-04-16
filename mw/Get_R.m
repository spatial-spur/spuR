function [R,DS] = Get_R(sig,qmax)
% Construct eigenvectors correpsonding to largest qmax eigenvalues of the
% nxn matrix sigma.
% Normalize the eigenvectors so that R'*R= I

n = size(sig,1);
[V,D] = eig(sig,'vector');
[D,ii] = sort(D,'descend');
DS = D(1:qmax);
R = V(:,ii(1:qmax));
% % Normalize so that R'*R/n = 1;
% rr = sum(R(:,1).^2);
% scl = sqrt(n/rr);
% R = scl*R;

end