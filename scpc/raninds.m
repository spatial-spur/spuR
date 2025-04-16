function [v] = raninds(n,M)
% yields capM+1 vector of random indices in {1,...,n} such that each pair is distinct
v = NaN(M+1,1);
j = floor(n*nextU);
for i = 1:M+1;
    v(i) = j+1;
    j=mod(j+1+floor(nextU()*(n-1)),n);
end;
    
end
