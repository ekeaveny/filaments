function [out] = apply_inverse_jacobian(J0invERROR_VECk,Cmat,Dmat,...
                                                               ERROR_VEC,k)
% APPLY_INVERSE_JACOBIAN   Calculates 
%                          [J_0^{-1} + \sum_{i=1}^{k} c_i d_i^T] f(X_k),
%                          i.e. Alg 2, Line 5 in the paper, as part of the
%                          bad Broyden's algorithm.
%
%

out = J0invERROR_VECk;

% Add on the  [\sum_{i=1}^{k} c_i d_i^T] f(X_k)  bit
for i = 1:k
    out = out + Cmat(:,i)*(Dmat(:,i)' * ERROR_VEC);
end

end

