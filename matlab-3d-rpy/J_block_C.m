function [C] = J_block_C(Filament,fac)
% J_BLOCK_C(Filament,fac)  Provides the terms related to the derivative of
%                          the velocity constraints with respect to the Lie
%                          algebra elements.

N_w = Filament.N_w;

Q = Filament.Q;

C = zeros(3*(N_w-1),3*N_w);

C(1:3,1:6) = fac*[rcross(QuaternionRotation(Q(1,1:4),[1;0;0])),...
                  rcross(QuaternionRotation(Q(2,1:4),[1;0;0]))];

for k=3:N_w
    
    i = 3*(k-2);
    
    C(i+1:i+3,1:i+3) = C(i-2:i,1:i+3);
    
    C(i+1:i+3,i+1:i+6) = C(i+1:i+3,i+1:i+6) ...
                 + fac*[rcross(QuaternionRotation(Q(k-1,1:4),[1;0;0])), ...
                        rcross(QuaternionRotation(Q(k,1:4),[1;0;0]))];
end

end

