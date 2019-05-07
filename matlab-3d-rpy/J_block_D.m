function [D] = J_block_D(Filament,mu)
% J_BLOCK_D(Filament,mu)   Provides the terms relating the velocity
%                          constraints to the Lagrange multipliers.

N_w = Filament.N_w;

D = zeros(3*(N_w-1));

for n=N_w:-1:1
    
    fac(n) = 1/(6*pi*mu*Filament.R(n));
    
end

for n=1:N_w-2
    
    i = 3*(n-1);
    j = i+3;
    
    D(i+1,j+1) = -fac(n+1);
    D(i+2,j+2) = -fac(n+1);
    D(i+3,j+3) = -fac(n+1);
    
    D(i+1,1) = fac(1);
    D(i+2,2) = fac(1);
    D(i+3,3) = fac(1);
    
end

D(end-2,1) = fac(1);
D(end-1,2) = fac(1);
D(end,3) = fac(1);

for i=1:3*(N_w-1)
    
    n = 1 + ceil(i/3);
    
    D(i,i) = D(i,i) + fac(n);
    
end

end

