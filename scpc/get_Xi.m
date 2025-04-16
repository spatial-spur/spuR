function[Xi] = get_Xi(X,m_vec)
   % Put elements of X into Xi
    n = size(m_vec,1);
    N = size(X,1);
    Xi = zeros(N,n);
    for l = 1:n;
        xl = get_xl_X(X,l,m_vec);
        if l == 1
           il = 0;
        else
           il = sum(m_vec(1:l-1,1));
        end
        ml = m_vec(l);
        Xi(il+1:il+ml,l) = xl;
    end;
    
end