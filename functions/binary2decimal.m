function x = binary2decimal( s )
%BINARY2DECIMAL Summary of this function goes here
%   Detailed explanation goes here
[m,n] = size(s);

% Convert to numbers
v = s - '0'; 
twos = pow2(n-1:-1:0);
toAdd = v.*twos;
toAdd(toAdd == 0) = [];

x = sum(toAdd);
end

