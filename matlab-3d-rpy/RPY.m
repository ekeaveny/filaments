function [] = RPY(Filaments,mu)
% RPY  Solves the Stokes flow problem using the Rotne-Prager-Yamakawa
%      tensor.
%
%   RPY(Filaments,mu)
%   sets the velocities and angular velocities of the Filament objects
%   for given forces and torques, for spherical particles, i, of radius 
%   Filament.R(i) at positions Filament.X(i) in an unbounded Newtonian 
%   fluid of viscosity mu.
%
%   Details of RPY solver: Wajnryb et al., 2013 Journal of Fluid Mechanics,
%   "Generalization of the Rotne-Prager-Yamakawa mobility and shear
%   disturbance tensors".
%
%   These expressions correpond to equations (23)-(26) in the paper.


N = length(Filaments);

for i=1:N
    
    Xi = Filaments(i).X;
    
    Ri = Filaments(i).R;
    
    V = zeros(3,Filaments(i).N_w);
    
    Omega = zeros(3,Filaments(i).N_w);
    
    for j=1:N
        
        Xj = Filaments(j).X;
        
        Fj = Filaments(j).F(1:3,:);
        
        Tj = Filaments(j).F(4:6,:);
        
        Rj = Filaments(j).R;
        
        for m=1:Filaments(i).N_w
            
            Xm = Xi(:,m);
            
            am = Ri(m);
            
            for n=1:Filaments(j).N_w
                
                Xn = Xj(:,n);
                
                r = ((Xm(1) - Xn(1))*(Xm(1) - Xn(1)) + (Xm(2) - Xn(2))*(Xm(2) - Xn(2)) + (Xm(3) - Xn(3))*(Xm(3) - Xn(3)))^0.5 + 10^-16;
                
                rhat = (Xm - Xn)/r;
                
                F = Fj(:,n);
                
                T = Tj(:,n);
                
                an = Rj(n);
                
                amax = max(am,an);
                
                amin = min(am,an);
                
                if r > (am + an) % The particles don't overlap.
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    % We begin with the translation-translation mobility.
                    
                    a0 = 1/(8*pi*mu*r);
                    a1 = 1 + (am*am + an*an)/(3*r*r);
                    a2 = 1 - (am*am + an*an)/(r*r);
                    
                    dot = F(1)*rhat(1) + F(2)*rhat(2) + F(3)*rhat(3);
                    
                    V(1,m) = V(1,m) + a0*(a1*F(1) + a2*dot*rhat(1));
                    V(2,m) = V(2,m) + a0*(a1*F(2) + a2*dot*rhat(2));
                    V(3,m) = V(3,m) + a0*(a1*F(3) + a2*dot*rhat(3));
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    % Next, we look at the rotation-rotation mobility.
                    
                    a0 = 1/(16*pi*mu*(r^3));
                    
                    dot = T(1)*rhat(1) + T(2)*rhat(2) + T(3)*rhat(3);
                    
                    Omega(1,m) = Omega(1,m) + a0*(3*dot*rhat(1) - T(1));
                    Omega(2,m) = Omega(2,m) + a0*(3*dot*rhat(2) - T(2));
                    Omega(3,m) = Omega(3,m) + a0*(3*dot*rhat(3) - T(3));
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    % Finally, we examine the translation-rotation mobilities.
                    
                    a0 = 1/(8*pi*mu*r*r);
                    
                    V(1,m) = V(1,m) + a0*(T(2)*rhat(3) - T(3)*rhat(2));
                    V(2,m) = V(2,m) + a0*(T(3)*rhat(1) - T(1)*rhat(3));
                    V(3,m) = V(3,m) + a0*(T(1)*rhat(2) - T(2)*rhat(1));
                    
                    Omega(1,m) = Omega(1,m) + a0*(F(2)*rhat(3) - F(3)*rhat(2));
                    Omega(2,m) = Omega(2,m) + a0*(F(3)*rhat(1) - F(1)*rhat(3));
                    Omega(3,m) = Omega(3,m) + a0*(F(1)*rhat(2) - F(2)*rhat(1));
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                elseif r > (amax - amin) % i.e. amax - amin < rij <= ai + aj
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    % We begin with the translation-translation mobility.
                    
                    a0 = 1/(6*pi*mu*am*an);
                    a1 = (16*(am + an)*r^3 - ((am - an)^2 + 3*r*r)^2)/(32*r^3);
                    a2 = (3*((am - an)^2 - r*r)^2)/(32*r^3);
                    
                    dot = F(1)*rhat(1) + F(2)*rhat(2) + F(3)*rhat(3);
                    
                    V(1,m) = V(1,m) + a0*(a1*F(1) + a2*dot*rhat(1));
                    V(2,m) = V(2,m) + a0*(a1*F(2) + a2*dot*rhat(2));
                    V(3,m) = V(3,m) + a0*(a1*F(3) + a2*dot*rhat(3));
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    % Next, we look at the rotation-rotation mobility.
                    
                    a0 = 1/(8*pi*mu*(am*an)^3);
                    A = (5*r^6 - 27*r^4*(am*am + an*an) + 32*r^3*(am^3 + an^3)...
                        -9*r*r*(am*am - an*an)^2 - (am - an)^4*(am*am + 4*am*an + an*an))/(64*r^3);
                    B = 3*((am - an)^2 - r*r)^2*(am*am + 4*am*an + an*an - r*r)/(64*r^3);
                    
                    dot = T(1)*rhat(1) + T(2)*rhat(2) + T(3)*rhat(3);
                    
                    Omega(1,m) = Omega(1,m) + a0*(B*dot*rhat(1) + A*T(1));
                    Omega(2,m) = Omega(2,m) + a0*(B*dot*rhat(2) + A*T(2));
                    Omega(3,m) = Omega(3,m) + a0*(B*dot*rhat(3) + A*T(3));
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    % Finally, we examine the translation-rotation mobilities.
                    
                    a0 = (1/(16*pi*mu*an*am^3))*((an*an + 2*an*(am + r) - 3*(am - r)^2)*(am - an + r)^2)/(8*r*r);
                    
                    V(1,m) = V(1,m) + a0*(T(2)*rhat(3) - T(3)*rhat(2));
                    V(2,m) = V(2,m) + a0*(T(3)*rhat(1) - T(1)*rhat(3));
                    V(3,m) = V(3,m) + a0*(T(1)*rhat(2) - T(2)*rhat(1));
                    
                    Omega(1,m) = Omega(1,m) + a0*(F(2)*rhat(3) - F(3)*rhat(2));
                    Omega(2,m) = Omega(2,m) + a0*(F(3)*rhat(1) - F(1)*rhat(3));
                    Omega(3,m) = Omega(3,m) + a0*(F(1)*rhat(2) - F(2)*rhat(1));
                    
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    
                    
                else % i.e. 0 <= rij <= amax - amin
                    
                    V(1,m) = V(1,m) + F(1)/(6*pi*mu*amax);
                    V(2,m) = V(2,m) + F(2)/(6*pi*mu*amax);
                    V(3,m) = V(3,m) + F(3)/(6*pi*mu*amax);
                    
                    Omega(1,m) = Omega(1,m) + T(1)/(8*pi*mu*amax^3);
                    Omega(2,m) = Omega(2,m) + T(2)/(8*pi*mu*amax^3);
                    Omega(3,m) = Omega(3,m) + T(3)/(8*pi*mu*amax^3);
                    
                    tr_fac = heaviside(am - an)*r/(8*pi*mu*am^3);
                    
                    V(1,m) = V(1,m) + tr_fac*(T(2)*rhat(3) - T(3)*rhat(2));
                    V(2,m) = V(2,m) + tr_fac*(T(3)*rhat(1) - T(1)*rhat(3));
                    V(3,m) = V(3,m) + tr_fac*(T(1)*rhat(2) - T(2)*rhat(1));
                    
                    Omega(1,m) = Omega(1,m) + tr_fac*(F(2)*rhat(3) - F(3)*rhat(2));
                    Omega(2,m) = Omega(2,m) + tr_fac*(F(3)*rhat(1) - F(1)*rhat(3));
                    Omega(3,m) = Omega(3,m) + tr_fac*(F(1)*rhat(2) - F(2)*rhat(1));
                    
                end % End overlap check.
                
            end
            
        end
        
    end
    
    Filaments(i).V = [V;Omega];
    
end

end % End function.

