
n = 5;
[x, y] = meshgrid(1:n, 1:n);
coords = [x(:), y(:)];
dist_mat = squareform(pdist(coords));

% Test different rho values
rho_vals = [0.001, 0.01, 0.05, 0.1];
results = zeros(size(rho_vals));

function [cbar] = getcbar(rhobar,distmat)
    % get cbar from rhobar and distance matrix
    % use bisection
    c0 = 10;
    c1 = 10;
    vd = lvech(distmat);  
    
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
      if c1 > 10000
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
    
    
% Run getcbar for each rho value
for i = 1:length(rho_vals)
    results(i) = getcbar(rho_vals(i), dist_mat);
end

% Print results
disp('MATLAB getcbar results:');
T = table(rho_vals', results', 'VariableNames', {'rho', 'cbar'});
disp(T);