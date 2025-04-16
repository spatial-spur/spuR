function [cv]=findcv(omega_mat,alpha,q)
% Find critical value given q
%   
global GQxw

% Find Critical Value: 
% Note this is cv/sqrt(q) from paper. sqrt(q) is included at the bottom
%  

% Get critical value of a given value of c .. begin with c0
icheck = 1;
itest = 0;

while itest == 0

  omega_mat_check = squeeze(omega_mat(1:q+1,1:q+1,icheck));

  % Find cv0 with rejection frequency >= alpha
  cv0 = tinv(1-alpha/2,q)/sqrt(q);  % Note cv0 is correct for c = infinity in mean case
 % icv0 = 0;
 % while icv0 == 0
 %   p = rejprob(cv0,omega_mat_check);
 %   if p >= alpha
 %     icv0 = 1;
 %   else
 %     cv0 = cv0/2;
 %   end
 %  end
      
  cv1 = cv0;
  i1 = 0;
  while i1 == 0
     p = rejprob(cv1,omega_mat_check);
     if p > alpha
        cv0 = cv1;
        cv1 = cv1 + 1/sqrt(q);
     else
         i1 = 1;
     end
  end
    
  % Bisection determination of cv
  cvm = (cv0+cv1)/2;
  while (cv1-cv0) > 0.001/sqrt(q)
      cvm = (cv0+cv1)/2;
      p = rejprob(cvm,omega_mat_check);
      if p > alpha
        cv0 = cvm;
      else
        cv1 = cvm;
      end
  end 
   
  % Find maximum rejection probability
  [p,ii]= maxrejprob(cvm,q,omega_mat);
  
  if ii == icheck
      itest = 1;
  else
      icheck = ii;
  end
    
end  
  
  cv = cvm.*sqrt(q);
 

end

