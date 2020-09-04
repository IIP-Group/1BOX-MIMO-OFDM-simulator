function t = dec2twos(x, nbits)

% DEC2TWOS Convert decimal integer to binary string two's complement.
% 
% Usage: T = DEC2TWOS(X, NBITS)
% 
% Converts the signed decimal integer given by X (either a scalar, vector,
% or matrix) to the two's complement representation as a string. If X is a
% vector or matrix, T(i, :) is the representation for X(i) (the shape of X
% is not preserved).
% 
% Example:
%     >> dec2twos([23 3 -23 -3])
%     
%     ans =
%     
%     010111
%     000011
%     101001
%     111101
% 
% Inputs:
%   -X: decimal integers to convert to two's complement.
%   -NBITS: number of bits in the representation 


x = x(:);

t = repmat('0', numel(x), nbits); % Initialize output:  Case for x = 0
if any(x > 0)
    t(x > 0, :) = dec2bin(x(x > 0), nbits);           % Case for x > 0
end
if any(x < 0)
    t(x < 0, :) = dec2bin(2^nbits + x(x < 0), nbits); % Case for x < 0
end