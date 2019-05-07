function [q] = qexp(u)
% QEXP(u)  Quaternion exponential map, exp(u) for an element of the Lie 
%          algebra, u.

theta = norm(u) + 10^-16;

q = zeros(1,4);

q(1) = cos(0.5*theta);

q(2:4) = sin(0.5*theta)*[u(1),u(2),u(3)]/theta;

end

