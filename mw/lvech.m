function [v] =  lvech(S)
% Computes lower triangle (below main diagonal) and puts in vector v
n = size(S,1);
v = NaN;
for i = 1:n-1;
    v = [v;S(i+1:end,i)];
end;
v = v(2:end,1);

end

