function [s] = jumble_s(s,m)
% randomly jumbles first m rows of location matrix s

n = size(s,1);
for i = 1:m;
    j = floor(nextU*n)+1;
    sj = s(j,:);
    s(j,:) = s(i,:);
    s(i,:) = sj;
end;

end

