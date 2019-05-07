function [F] = J_block_F(Filament,dt,mu)
% J_BLOCK_F(Filament,dt,mu)   Provides the terms relating the Lie algebra 
%                             update equations to the Lagrange multipliers.

N_w = Filament.N_w;

DL = Filament.DL;

R = Filament.R;

U = Filament.U;

Q = Filament.Q;

F = zeros(3*N_w,3*(N_w-1));

for i=1:N_w-1
    
    pos = 3*(i-1);
    
    fac = dt*DL/(24*pi*mu*R(i)^3);
    
    tmat = rcross(QuaternionRotation(Q(i,1:4),[1;0;0]))';
    umat = rcross(U(1:3,i))';
    
    F(pos+1:pos+3,pos+1:pos+3) = fac*(tmat - 0.5*umat*tmat + umat^2*tmat/12);
    
    fac = dt*DL/(24*pi*mu*R(i+1)^3);
    
    tmat = rcross(QuaternionRotation(Q(i+1,1:4),[1;0;0]))';
    umat = rcross(U(1:3,i+1))';
    
    F(pos+4:pos+6,pos+1:pos+3) = fac*(tmat - 0.5*umat*tmat + umat^2*tmat/12);
    
end

end

