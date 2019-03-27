function [VX,VY,VZ,OMEGAX,OMEGAY,OMEGAZ] = RPY(FX,FY,FZ,TAUX,TAUY,TAUZ,X,Y,Z,a,mu)
% RPY  Solves the Stokes flow problem using the Rotne-Prager-Yamakawa
%      tensor.
%
%   [VX,VY,VZ,OMEGAX,OMEGAY,OMEGAZ] = RPY(FX,FY,FZ, ...
%                                             TAUX,TAUY,TAUZ,X,Y,Z,a,mu)
%   returns velocities and angular velocities for given forces and torques
%   for spherical particles of radius a at positions [X,Y,Z] in an 
%   unbounded Newtonian fluid of viscosity mu.
%
%   Details of RPY solver: Wajnryb et al., 2013 Journal of Fluid Mechanics,
%   "Generalization of the Rotne-Prager-Yamakawa mobility and shear
%   disturbance tensors".
%
%   These expressions correpond to equations (23)-(26) in the paper.

N = size(FX,1); % Extract the system size.

% Initialise the output arrays using the self mobilities:

FC = 6*pi*a*mu;
TC = 8*pi*mu*(a^3);

VX = FX/FC; VY = FY/FC; VZ = FZ/FC;
OMEGAX = TAUX/TC; OMEGAY = TAUY/TC; OMEGAZ = TAUZ/TC;

% Loop over the particles and evaluate the remaining velocities:

for jj=1:N
    for kk=1+jj:N
        
        % Define required quantities:
        
        rij = norm([X(jj),Y(jj),Z(jj)] - [X(kk),Y(kk),Z(kk)]);
        
        rhat = ([X(jj),Y(jj),Z(jj)] - [X(kk),Y(kk),Z(kk)])/rij;
        
        % We implement different expressions for the mobilities depending
        % on whether the particles overlap or not:
        
        if rij > (2*a) % The particles don't overlap.
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            % We begin with the translation-translation mobility.
            
            a0 = 1/(8*pi*mu*rij);
            a1 = 1 + (2*a*a)/(3*rij*rij);
            a2 = 1 - (2*a*a)/(rij*rij);
            
            % Effect of particle kk on jj.
            
            rhat_dot_force = rhat(1)*FX(kk) + rhat(2)*FY(kk) + rhat(3)*FZ(kk);
            
            VX(jj) = VX(jj) + a0*(a1*FX(kk) + a2*rhat_dot_force*rhat(1));
            VY(jj) = VY(jj) + a0*(a1*FY(kk) + a2*rhat_dot_force*rhat(2));
            VZ(jj) = VZ(jj) + a0*(a1*FZ(kk) + a2*rhat_dot_force*rhat(3));
            
            % Effect of particle jj on kk.
            
            rhat_dot_force = rhat(1)*FX(jj) + rhat(2)*FY(jj) + rhat(3)*FZ(jj);
            
            VX(kk) = VX(kk) + a0*(a1*FX(jj) + a2*rhat_dot_force*rhat(1));
            VY(kk) = VY(kk) + a0*(a1*FY(jj) + a2*rhat_dot_force*rhat(2));
            VZ(kk) = VZ(kk) + a0*(a1*FZ(jj) + a2*rhat_dot_force*rhat(3));
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            % Next, we look at the rotation-rotation mobility.
            
            a0 = 1/(16*pi*mu*(rij^3));
            
            % Effect of particle kk on jj.
            
            rhat_dot_torque = rhat(1)*TAUX(kk) + rhat(2)*TAUY(kk) + rhat(3)*TAUZ(kk);
            
            OMEGAX(jj) = OMEGAX(jj) + a0*(3*rhat_dot_torque*rhat(1) - TAUX(kk));
            OMEGAY(jj) = OMEGAY(jj) + a0*(3*rhat_dot_torque*rhat(2) - TAUY(kk));
            OMEGAZ(jj) = OMEGAZ(jj) + a0*(3*rhat_dot_torque*rhat(3) - TAUZ(kk));
            
            % Effect of particle jj on kk.
            
            rhat_dot_torque = rhat(1)*TAUX(jj) + rhat(2)*TAUY(jj) + rhat(3)*TAUZ(jj);
            
            OMEGAX(kk) = OMEGAX(kk) + a0*(3*rhat_dot_torque*rhat(1) - TAUX(jj));
            OMEGAY(kk) = OMEGAY(kk) + a0*(3*rhat_dot_torque*rhat(2) - TAUY(jj));
            OMEGAZ(kk) = OMEGAZ(kk) + a0*(3*rhat_dot_torque*rhat(3) - TAUZ(jj));
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            % Finally, we examine the translation-rotation mobilities.
            
            a0 = 1/(8*pi*mu*rij*rij);
            
            % Effect of particle kk on jj.
            
            VX(jj) = VX(jj) + a0*(TAUY(kk)*rhat(3) - TAUZ(kk)*rhat(2));
            VY(jj) = VY(jj) + a0*(TAUZ(kk)*rhat(1) - TAUX(kk)*rhat(3));
            VZ(jj) = VZ(jj) + a0*(TAUX(kk)*rhat(2) - TAUY(kk)*rhat(1));
            
            OMEGAX(jj) = OMEGAX(jj) + a0*(FY(kk)*rhat(3) - FZ(kk)*rhat(2));
            OMEGAY(jj) = OMEGAY(jj) + a0*(FZ(kk)*rhat(1) - FX(kk)*rhat(3));
            OMEGAZ(jj) = OMEGAZ(jj) + a0*(FX(kk)*rhat(2) - FY(kk)*rhat(1));
            
            % Effect of particle jj on kk. Note that for this term we have
            % a dependence on the sign of the vector rhat.
            
            VX(kk) = VX(kk) - a0*(TAUY(jj)*rhat(3) - TAUZ(jj)*rhat(2));
            VY(kk) = VY(kk) - a0*(TAUZ(jj)*rhat(1) - TAUX(jj)*rhat(3));
            VZ(kk) = VZ(kk) - a0*(TAUX(jj)*rhat(2) - TAUY(jj)*rhat(1));
            
            OMEGAX(kk) = OMEGAX(kk) - a0*(FY(jj)*rhat(3) - FZ(jj)*rhat(2));
            OMEGAY(kk) = OMEGAY(kk) - a0*(FZ(jj)*rhat(1) - FX(jj)*rhat(3));
            OMEGAZ(kk) = OMEGAZ(kk) - a0*(FX(jj)*rhat(2) - FY(jj)*rhat(1));
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        else % The particles overlap.
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            % We begin with the translation-translation mobility.
            
            a0 = 1/FC;
            a1 = 1 - (9*rij)/(32*a);
            a2 = (3*rij)/(32*a);
            
            % Effect of particle kk on jj.
            
            rhat_dot_force = rhat(1)*FX(kk) + rhat(2)*FY(kk) + rhat(3)*FZ(kk);
            
            VX(jj) = VX(jj) + a0*(a1*FX(kk) + a2*rhat_dot_force*rhat(1));
            VY(jj) = VY(jj) + a0*(a1*FY(kk) + a2*rhat_dot_force*rhat(2));
            VZ(jj) = VZ(jj) + a0*(a1*FZ(kk) + a2*rhat_dot_force*rhat(3));
            
            % Effect of particle jj on kk.
            
            rhat_dot_force = rhat(1)*FX(jj) + rhat(2)*FY(jj) + rhat(3)*FZ(jj);
            
            VX(kk) = VX(kk) + a0*(a1*FX(jj) + a2*rhat_dot_force*rhat(1));
            VY(kk) = VY(kk) + a0*(a1*FY(jj) + a2*rhat_dot_force*rhat(2));
            VZ(kk) = VZ(kk) + a0*(a1*FZ(jj) + a2*rhat_dot_force*rhat(3));
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            % Next, we look at the rotation-rotation mobility.
            
            a0 = 1/TC;
            a1 = (9*rij)/(32*a) - (3/64)*((rij/a)^3);
            a2 = 1 - (27*rij)/(32*a) + (5/64)*((rij/a)^3);
            
            % Effect of particle kk on jj.
            
            rhat_dot_torque = rhat(1)*TAUX(kk) + rhat(2)*TAUY(kk) + rhat(3)*TAUZ(kk);
            
            OMEGAX(jj) = OMEGAX(jj) + a0*(a1*rhat_dot_torque*rhat(1) + a2*TAUX(kk));
            OMEGAY(jj) = OMEGAY(jj) + a0*(a1*rhat_dot_torque*rhat(2) + a2*TAUY(kk));
            OMEGAZ(jj) = OMEGAZ(jj) + a0*(a1*rhat_dot_torque*rhat(3) + a2*TAUZ(kk));
            
            % Effect of particle jj on kk.
            
            rhat_dot_torque = rhat(1)*TAUX(jj) + rhat(2)*TAUY(jj) + rhat(3)*TAUZ(jj);
            
            OMEGAX(kk) = OMEGAX(kk) + a0*(a1*rhat_dot_torque*rhat(1) + a2*TAUX(jj));
            OMEGAY(kk) = OMEGAY(kk) + a0*(a1*rhat_dot_torque*rhat(2) + a2*TAUY(jj));
            OMEGAZ(kk) = OMEGAZ(kk) + a0*(a1*rhat_dot_torque*rhat(3) + a2*TAUZ(jj));
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            % Finally, we examine the translation-rotation mobilities.
            
            a0 = (1/(16*pi*mu*a*a))*((rij/a) - (3*rij*rij)/(8*a*a));
            
            % Effect of particle kk on jj.
            
            VX(jj) = VX(jj) + a0*(TAUY(kk)*rhat(3) - TAUZ(kk)*rhat(2));
            VY(jj) = VY(jj) + a0*(TAUZ(kk)*rhat(1) - TAUX(kk)*rhat(3));
            VZ(jj) = VZ(jj) + a0*(TAUX(kk)*rhat(2) - TAUY(kk)*rhat(1));
            
            OMEGAX(jj) = OMEGAX(jj) + a0*(FY(kk)*rhat(3) - FZ(kk)*rhat(2));
            OMEGAY(jj) = OMEGAY(jj) + a0*(FZ(kk)*rhat(1) - FX(kk)*rhat(3));
            OMEGAZ(jj) = OMEGAZ(jj) + a0*(FX(kk)*rhat(2) - FY(kk)*rhat(1));
            
            % Effect of particle jj on kk. Note that for this term we have
            % a dependence on the sign of the vector rhat.
            
            VX(kk) = VX(kk) - a0*(TAUY(jj)*rhat(3) - TAUZ(jj)*rhat(2));
            VY(kk) = VY(kk) - a0*(TAUZ(jj)*rhat(1) - TAUX(jj)*rhat(3));
            VZ(kk) = VZ(kk) - a0*(TAUX(jj)*rhat(2) - TAUY(jj)*rhat(1));
            
            OMEGAX(kk) = OMEGAX(kk) - a0*(FY(jj)*rhat(3) - FZ(jj)*rhat(2));
            OMEGAY(kk) = OMEGAY(kk) - a0*(FZ(jj)*rhat(1) - FX(jj)*rhat(3));
            OMEGAZ(kk) = OMEGAZ(kk) - a0*(FX(jj)*rhat(2) - FY(jj)*rhat(1));
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        end % End overlap check.
        
    end % End inner particle loop.
end % End outer particle loop.

end % End function.

