function [p] = GQp(eta,GQxw)
%  Compute probability by Gaussian quadrature
%  GQxw .. ordinates and weights
%  eta .. eta values from formula
u = GQxw(:,1);
w = GQxw(:,2);
a = 1-u.^2;
b = 1+repmat(eta',size(u,1),1)./repmat(a,1,size(eta,1));
c = 1./sqrt(a.*prod(b,2));
p = (2/pi)*(w'*c);

end

