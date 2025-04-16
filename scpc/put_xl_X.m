function [X] = put_xl_X(X,xl,l,m_vec)
    % put submatrix xl into X
    ml = m_vec(l);
    if l == 1
       X(1:ml,:) = xl;
    else
       ii = sum(m_vec(1:l-1,1));
       X(ii+1:ii+ml,:) = xl;
    end
   
end
  
