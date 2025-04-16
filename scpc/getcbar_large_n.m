function [cbar] = getcbar_large_n(rhobar,distmat_mN)
% get cbar from rhobar and distance matrix
% use bisection
c0 = 10;
c1 = 10;
N = size(distmat_mN,3);
vd = lvech(squeeze(distmat_mN(:,:,1)));
for i = 2:N;
    vd = [vd;lvech(squeeze(distmat_mN(:,:,i)))];
end;

% Adjust so c0 yield v > rhobar and c1 yields v < rhobar
i1 = 0;
jj = 0;
while i1 == 0
  v = mean(exp(-c0*vd));
  i1 = v > rhobar;
  if i1 == 0
      c1 = c0;
      c0 = c0/2;
      jj = jj+1;
  end
  if jj > 500
   error('rhobar too large');
  end
end

% Verify that c1 yields value lt rhobar
i1 = 0;
jj = 0;
while i1 == 0
  v = mean(exp(-c1*vd));
  i1 = v < rhobar;
  if i1 == 0
      c0 = c1;
      c1 = 2*c1;
      jj = jj+1;
  end
  if c1 > 5000
      i1 = 1;
  end
  if jj > 500
   error('rhobar too small');
  end
end

% Bisection determination of cbar
while (c1-c0) > 0.001
   cm = sqrt(c0*c1);
   v = mean(exp(-cm*vd));
   if v < rhobar
      c1 = cm;
   else
      c0 = cm;
   end
end

cbar = sqrt(c0*c1);


end

