function g = get_ha_parm_I0(om_ho,om_i0,om_bm,e)
% Find value of parameter that gives (roughly) 50% power
pow = 1;
gtry = 1;
while pow > 0.5
  g = gtry;
  % Verify that power is less than 50%;
  pow = getpow_qf(om_ho,om_i0+g*om_bm,e);
  % [g pow]
  gtry = g/2;
end
g1 = g;

pow = 0;
gtry = 30;
while pow < 0.5
  g = gtry;
  % Verify that power is less than 50%;
  pow = getpow_qf(om_ho,om_i0+g*om_bm,e);
  % [pow g]
  gtry = g*2;
end
g2 = g;

ii = 1;
while (abs(pow-0.5) > 0.01)
    g = (g1+g2)/2;
    pow = getpow_qf(om_ho,om_i0+g*om_bm,e);
    if pow > 0.5
        g2 = g;
    elseif pow < 0.5
        g1 = g;
    end
    % [pow g1 g2]
    ii = ii+1;
    if (ii > 20) 
        break 
    end
end

end