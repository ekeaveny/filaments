function [E] = J_block_E(Filament,dt,mu)
% J_BLOCK_E(Filament,dt,mu)   Provides the terms relating the Lie algebra 
%                             update equations on the Lie algebra elements.


N_w = Filament.N_w;

E = zeros(3*N_w);

DL = Filament.DL;

KB = Filament.KB;

KT = Filament.KT;

strain_twist = Filament.StrainTwist;

Q = Filament.Q;

Lambda = Filament.Lambda;

U = Filament.U;

Omega = Filament.V(4:6,:);

R = Filament.R;

% We deal with the end-particle cases seperately.

% n = 1

Tfac = 1/(8*pi*mu*R(1)^3);

t = QuaternionRotation(Q(1,1:4),[1;0;0]);
n = QuaternionRotation(Q(1,1:4),[0;1;0]);
b = QuaternionRotation(Q(1,1:4),[0;0;1]);

tp1 = QuaternionRotation(Q(2,1:4),[1;0;0]);
np1 = QuaternionRotation(Q(2,1:4),[0;1;0]);
bp1 = QuaternionRotation(Q(2,1:4),[0;0;1]);

% Diagonal block

constraint_part = -0.5*DL*rcross(Lambda(:,1))*rcross(t);

beta = -strain_twist(3) + 0.5*(np1' * b - n' * bp1)/DL;

elastic_part = 0.5*(t + tp1)*(cross(b,np1) - cross(n,bp1))'/DL;

elastic_part = 0.5*KT*(beta*rcross(t) + elastic_part);

elastic_part = elastic_part + KB*(rcross(tp1)*rcross(t)/DL...
    - 0.5*strain_twist(1)*rcross(n)...
    - 0.5*strain_twist(2)*rcross(b));

dOmega = Tfac*(elastic_part + constraint_part);

Lie_mat = rcross(U(1:3,1))';
Omega_mat = rcross(Omega(:,1));

block = -2*dt*(dOmega - 0.5*(Omega_mat + Lie_mat*dOmega) + ...
    (rcross(cross(U(1:3,1),Omega(:,1))) + Lie_mat*Omega_mat + Lie_mat^2*dOmega)/12)/3;

E(1:3,1:3) = block;

% Off-diagonal block

elastic_part = 0.5*(t + tp1)*(cross(b,np1) - cross(n,bp1))'/DL;

elastic_part = 0.5*KT*(beta*rcross(tp1) - elastic_part);

elastic_part = elastic_part + KB*(rcross(t)'*rcross(tp1)/DL...
    - 0.5*strain_twist(1)*rcross(np1)...
    - 0.5*strain_twist(2)*rcross(bp1));

dOmega = Tfac*elastic_part;

block = -2*dt*(dOmega - 0.5*Lie_mat*dOmega + Lie_mat^2*dOmega/12)/3;

E(1:3,4:6) = block;

% n = N_w

Tfac = 1/(8*pi*mu*R(N_w)^3);

tm1 = QuaternionRotation(Q(N_w-1,1:4),[1;0;0]);
nm1 = QuaternionRotation(Q(N_w-1,1:4),[0;1;0]);
bm1 = QuaternionRotation(Q(N_w-1,1:4),[0;0;1]);

t = QuaternionRotation(Q(N_w,1:4),[1;0;0]);
n = QuaternionRotation(Q(N_w,1:4),[0;1;0]);
b = QuaternionRotation(Q(N_w,1:4),[0;0;1]);

% Diagonal block

constraint_part = -0.5*DL*rcross(Lambda(:,end))*rcross(t);

beta = -strain_twist(3) + 0.5*(n' * bm1 - nm1' * b)/DL;

elastic_part = 0.5*(tm1 + t)*(cross(bm1,n) - cross(nm1,b))'/DL;

elastic_part = -0.5*KT*(beta*rcross(t) - elastic_part);

elastic_part = elastic_part - KB*(rcross(tm1)'*rcross(t)/DL...
    - 0.5*strain_twist(1)*rcross(n)...
    - 0.5*strain_twist(2)*rcross(b));

dOmega = Tfac*(elastic_part + constraint_part);

Lie_mat = rcross(U(1:3,N_w))';
Omega_mat = rcross(Omega(:,N_w));

block = -2*dt*(dOmega - 0.5*(Omega_mat + Lie_mat*dOmega) + ...
    (rcross(cross(U(1:3,N_w),Omega(:,N_w))) + Lie_mat*Omega_mat + Lie_mat^2*dOmega)/12)/3;

E(end-2:end,end-2:end) = block;

% Off-diagonal block

elastic_part = 0.5*(tm1 + t)*(cross(bm1,n) - cross(nm1,b))'/DL;

elastic_part = -0.5*KT*(beta*rcross(tm1) + elastic_part);

elastic_part = elastic_part - KB*(rcross(t)*rcross(tm1)/DL...
    - 0.5*strain_twist(1)*rcross(nm1)...
    - 0.5*strain_twist(2)*rcross(bm1));

dOmega = Tfac*elastic_part;

block = -2*dt*(dOmega - 0.5*Lie_mat*dOmega + Lie_mat^2*dOmega/12)/3;

E(end-2:end,end-5:end-3) = block;

% Now for the interior points, which experience an elastic moment and a
% constraint torque on each side, giving their entries the most general
% form.

for k=2:N_w-1
    
    Tfac = 1/(8*pi*mu*R(k)^3);
    
    tm1 = QuaternionRotation(Q(k-1,1:4),[1;0;0]);
    nm1 = QuaternionRotation(Q(k-1,1:4),[0;1;0]);
    bm1 = QuaternionRotation(Q(k-1,1:4),[0;0;1]);
    
    t = QuaternionRotation(Q(k,1:4),[1;0;0]);
    n = QuaternionRotation(Q(k,1:4),[0;1;0]);
    b = QuaternionRotation(Q(k,1:4),[0;0;1]);
    
    tp1 = QuaternionRotation(Q(k+1,1:4),[1;0;0]);
    np1 = QuaternionRotation(Q(k+1,1:4),[0;1;0]);
    bp1 = QuaternionRotation(Q(k+1,1:4),[0;0;1]);
    
    Lie_mat = rcross(U(1:3,k))';
    Omega_mat = rcross(Omega(:,k));
    
    % Construct the diagonal block first.
    
    constraint_part = -0.5*DL*rcross(Lambda(:,k-1) + Lambda(:,k))*rcross(t);
    
    beta = -strain_twist(3) + 0.5*(np1' * b...
        - n' * bp1)/DL;
    
    elastic_part_right = 0.5*(t + tp1)...
        *(cross(b,np1)...
        - cross(n,bp1))'/DL;
    
    elastic_part_right = 0.5*KT*(beta*rcross(t) + elastic_part_right);
    
    elastic_part_right = elastic_part_right + KB*(rcross(tp1)*rcross(t)/DL...
        - 0.5*strain_twist(1)*rcross(n)...
        - 0.5*strain_twist(2)*rcross(b));
    
    beta = -strain_twist(3) + 0.5*(n' * bm1...
        - nm1' * b)/DL;
    
    elastic_part_left = 0.5*(tm1 + t)...
        *(cross(bm1,n)...
        - cross(nm1,b))'/DL;
    
    elastic_part_left = -0.5*KT*(beta*rcross(t) - elastic_part_left);
    
    elastic_part_left = elastic_part_left - KB*(rcross(tm1)'*rcross(t)/DL...
        - 0.5*strain_twist(1)*rcross(n)...
        - 0.5*strain_twist(2)*rcross(b));
    
    dOmega = Tfac*(constraint_part + elastic_part_right + elastic_part_left);
    
    block = -2*dt*(dOmega - 0.5*(Omega_mat + Lie_mat*dOmega) + ...
        (rcross(cross(U(1:3,k),Omega(:,k))) + Lie_mat*Omega_mat + Lie_mat^2*dOmega)/12)/3;
    
    E(3*(k-1)+1:3*(k-1)+3,3*(k-1)+1:3*(k-1)+3) = block;
    
    % Now for the off-diagonal blocks, starting with the left one.
    
    elastic_part = 0.5*(tm1 + t)...
        *(cross(bm1,n)...
        - cross(nm1,b))'/DL;
    
    elastic_part = -0.5*KT*(beta*rcross(tm1) + elastic_part);
    
    elastic_part = elastic_part - KB*(rcross(t)*rcross(tm1)/DL...
        - 0.5*strain_twist(1)*rcross(nm1)...
        - 0.5*strain_twist(2)*rcross(bm1));
    
    dOmega = Tfac*elastic_part;
    
    block = -2*dt*(dOmega - 0.5*Lie_mat*dOmega + Lie_mat^2*dOmega/12)/3;
    
    E(3*(k-1)+1:3*(k-1)+3,3*(k-2)+1:3*(k-2)+3) = block;
    
    % Finally, the right off-diagonal block.
    
    elastic_part = 0.5*(t + tp1)...
        *(cross(b,np1)...
        - cross(n,bp1))'/DL;
    
    elastic_part = 0.5*KT*(beta*rcross(tp1) - elastic_part);
    
    elastic_part = elastic_part + KB*(rcross(t)'*rcross(tp1)/DL...
        - 0.5*strain_twist(1)*rcross(np1)...
        - 0.5*strain_twist(2)*rcross(bp1));
    
    dOmega = Tfac*elastic_part;
    
    block = -2*dt*(dOmega - 0.5*Lie_mat*dOmega + Lie_mat^2*dOmega/12)/3;
    
    E(3*(k-1)+1:3*(k-1)+3,3*(k)+1:3*(k)+3) = block;
    
end

% All that remains is to add on the identity part.

for i=1:3*N_w
    
    E(i,i) = E(i,i) + 1;
    
end

end

