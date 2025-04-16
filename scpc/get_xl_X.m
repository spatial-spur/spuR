function [xl] = get_xl_X(X,l,m_vec)
    % get submatrix xl from X
    ml = m_vec(l);
    if l == 1
       xl=X(1:ml,:);
    else
       ii = sum(m_vec(1:l-1,1));
       xl = X(ii+1:ii+ml,:);
    end
   
end