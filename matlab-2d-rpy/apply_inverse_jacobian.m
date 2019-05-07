function [DeltaX] = apply_inverse_jacobian(J0invERROR_VEC,Cmat,Dmat,...
                                                               ERROR_VEC,k)
% APPLY_INVERSE_JACOBIAN   Calculates 
%                          [J_0^{-1} + \sum_{i=1}^{k} c_i d_i^T] f(X_k),
%                          i.e. Alg 2, Line 5 in the paper, as part of the
%                          bad Broyden's algorithm.
%
% [DeltaX] = apply_inverse_jacobian(J0invERROR_VEC,C_Mat,D_Mat,ERROR_VEC,k)

DeltaX = J0invERROR_VEC;

for i = 1:k
    DeltaX = DeltaX + Cmat(:,i)*(Dmat(:,i)'*ERROR_VEC);
end

end