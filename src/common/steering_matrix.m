function A = steering_matrix(p, thetaDeg)
%STEERING_MATRIX Sparse-array steering matrix.
% Convention used everywhere in this package:
%   a_m(theta) = exp(-j*pi*p_m*sin(theta))
% where p is in units of lambda/2.

p = p(:);
thetaDeg = thetaDeg(:).';
A = exp(-1j*pi*(p * sind(thetaDeg)));
end
