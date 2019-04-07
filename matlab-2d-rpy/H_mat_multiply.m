function [DeltaX] = H_mat_multiply(J0invERROR_VEC,C_Mat,D_Mat,ERROR_VEC,k)
% H_MAT_MULTIPLY   Calculates [J_0^{-1} + \sum_{i=1}^{k} c_i d_i^T] f(X_k),
%                  i.e. Alg 2, Line 5 in the paper, as part of the bad
%                  Broyden's algorithm.
%
%    [DeltaX] = H_mat_multiply(J0invERROR_VEC,C_Mat,D_Mat,ERROR_VEC,k)

DeltaX = J0invERROR_VEC;

for i = 1:k
    DeltaX = DeltaX + C_Mat(:,i)*(D_Mat(:,i)'*ERROR_VEC);
end

end