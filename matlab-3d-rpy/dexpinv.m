function [out] = dexpinv(u,v)
% DEXPINV  dexpinv equation for (R^3, \times).
%          Equation (58) in paper, where v = Omega.
%
%   dexpinv(u,v) = dexp_u^{-1}(v)

theta = norm(u);

if theta < 10^-6
    fac = -1/12;
else
    fac = (0.5*theta*cot(0.5*theta) - 1)/(theta^2);
end

out = v - 0.5*cross(u,v) - fac*cross(u,cross(u,v));

end

