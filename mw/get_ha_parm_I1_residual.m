function c = get_ha_parm_I1_residual(om_ho,distmat,R,e,M)

% Find value of parameter that gives (roughly) 50% power
pow = 1;
ctry = getcbar(0.95,distmat);
while pow > 0.5
  c = ctry;
  sigdm_c = get_sigma_residual(distmat,c,M);
  om_c = R'*sigdm_c*R;
  pow = getpow_qf(om_ho,om_c,e);
  % Verify that power is less than 50%;
  %[ctry pow]
  ctry = ctry/2;
end
c1 = c;
om_c1 = om_c;

pow = 0;
ctry = getcbar(0.01,distmat);
while pow < 0.5
  c = ctry;
  sigdm_c = get_sigma_residual(distmat,c,M);
  om_c = R'*sigdm_c*R;
  pow = getpow_qf(om_ho,om_c,e);
  % Verify that power is less than 50%;
  %[ctry pow]
  ctry = 2*ctry;
end
c2 = c;
om_c2 = om_c;

ii = 1;
while (abs(pow-0.5) > 0.01)
    c = (c1+c2)/2;
    sigdm_c = get_sigma_residual(distmat,c,M);
    om_c = R'*sigdm_c*R;
    pow = getpow_qf(om_ho,om_c,e);
    if pow > 0.5
        c2 = c;
    elseif pow < 0.5
        c1 = c;
    end
    %[pow c1 c2]
    ii = ii+1;
    if (ii > 20) 
        break 
    end
end

end