function [u] = nextU()
% Random number generator
global random_t

random_t=mod(64389*random_t + 1,2^32);
u = random_t/(2^32);

end

