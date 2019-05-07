function [product] = QuaternionProduct(p,q)
% QUATERNIONPRODUCT(p,q)   Returns p.q for two quaternions p and q.

product = zeros(1,4);

product(1) = p(1)*q(1) - p(2)*q(2) - p(3)*q(3) - p(4)*q(4);
product(2) = p(1)*q(2) + q(1)*p(2) + p(3)*q(4) - p(4)*q(3);
product(3) = p(1)*q(3) + q(1)*p(3) + p(4)*q(2) - p(2)*q(4);
product(4) = p(1)*q(4) + q(1)*p(4) + p(2)*q(3) - p(3)*q(2);

end

