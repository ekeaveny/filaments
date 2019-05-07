function [x] = blockwise_backslash_jacobian(Filaments,ERROR_VEC)
% BLOCKWISE_BACKSLASH_JACOBIAN   Calculates  [J_0^{-1}] f(X_k).
%
%   x = BLOCKWISE_BACKSLASH_JACOBIAN(Filaments,b) returns x = J_0\b where
%   J_0 is the approximate Jacobian of the system. J_0 is stored blockwise
%   in the Filament object, with the full approximate Jacobian being of the
%   form J0 = [A1  O  O ..
%               O A2  O ..
%               O  O A3 ..
%               :  :  :   ]. 
%
%   Each block A1, A2, etc., is the Jacobian of the system for each indivi-
%   dual filament. Solving x = [J_0^{-1}] f(X_k)  is equivalent, therefore, 
%   to solving multiple filament-sized systems.
%
for i = length(Filaments):-1:1
    N(i) = 6*Filaments(i).N_w;
end
N = [0,cumsum(N)];
x = zeros(N(end),1);

% Calculate [J_0^{-1}] f(X_k) blockwise
% We have already done the LU-decomp on the filament-sized J_0 blocks,
% so this is easy.
for i = 1:length(Filaments)
    x(N(i)+1:N(i+1)) = Filaments(i).InvertLocalBlock(...
                                                ERROR_VEC(N(i)+1:N(i+1)) );
end

end

