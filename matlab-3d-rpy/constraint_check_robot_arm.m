function [concheck,ERROR_VEC] = constraint_check_robot_arm(...
                                                       Filaments,dt,nt,tol)
% CONSTRAINT_CHECK_ROBOT_ARM  Creates the error vector for the Broyden's
%                             step and finds whether it is less than a 
%                             specified tolerance.
%
%   This corresponds to the system f(X*) at equation (62) in the paper.
%
%   [concheck, ERROR_VEC] = constraint_check_robot_arm(...
%                                                      Filaments,dt,nt,tol)
%   takes in, through Filament objects,
%      the calculated next timestep positions/orientations (Filaments.X),
%                  current timestep positions/orientations (Filaments.Xm1),
%             and previous timestep positions/orientations (Filaments.Xm2),
%             the calculated velocities/angular velocities (Filaments.V);
%   timestep dt, timestep number nt, 
%   and tolerance (epsilon in paper) tol.
%                                                                  
%   It returns ERROR_VEC (f(X) in the paper) and whether ||f(X)|| < tol.
%   concheck = 0 means ||f(X)|| < tol. Otherwise, concheck = 1.

concheck = 0;

N_sw = length(Filaments);

c_1 = 4.0/3.0;
c_2 = 1.0/3.0;
c_4 = 2.0/3.0;

for j_w=N_sw:-1:1
    N(1,j_w) = 6*Filaments(j_w).N_w;
end

N = [0,cumsum(N)];
Np = N(end);

ERROR_VEC = zeros(Np,1);

for n=1:N_sw
    
    V = Filaments(n).V;
    
    X   = Filaments(n).X;   % Positions of filament segments at next t'step
    Xm1 = Filaments(n).Xm1; % ... at current timestep
    Xm2 = Filaments(n).Xm2; % ... at previous timestep
    
    U = Filaments(n).U;     % Lie algebra elements associated with the 
                            % orientation of the constituent particles.
    
    if(nt == 1) % Euler   
        for j_w = 1:Filaments(n).N_w % segment number
            
            ERROR1 = X(:,j_w) - Xm1(:,j_w) - dt*V(1:3,j_w);
            
            ERROR3 = U(1:3,j_w) - dt*dexpinv(U(1:3,j_w),V(4:6,j_w));

            loc1 = N(n) + 3*(j_w - 1 + (j_w ~= 1)*Filaments(n).N_w);
            ERROR_VEC(loc1+1) = ERROR1(1);
            ERROR_VEC(loc1+2) = ERROR1(2);
            ERROR_VEC(loc1+3) = ERROR1(3);
            
            loc3 = N(n) + 3*j_w;
            ERROR_VEC(loc3+1) = ERROR3(1);
            ERROR_VEC(loc3+2) = ERROR3(2);
            ERROR_VEC(loc3+3) = ERROR3(3);

            if(max(abs(ERROR1)) > tol || max(abs(ERROR3)) > tol)
                concheck = 1;
            end
        end
        
    else %% BDF2
        for j_w = 1:Filaments(n).N_w
            
            ERROR1 = X(:,j_w) - c_1*Xm1(:,j_w) ...
                           + c_2*Xm2(:,j_w) - dt*c_4*V(1:3,j_w);
                       
            ERROR3 = U(1:3,j_w) - c_2*U(4:6,j_w) ...
                     - c_4*dt*dexpinv(U(1:3,j_w),V(4:6,j_w));
            
            loc1 = N(n) + 3*(j_w - 1 + (j_w ~= 1)*Filaments(n).N_w);
            ERROR_VEC(loc1+1) = ERROR1(1);
            ERROR_VEC(loc1+2) = ERROR1(2);
            ERROR_VEC(loc1+3) = ERROR1(3);
            
            loc3 = N(n) + 3*j_w;
            ERROR_VEC(loc3+1) = ERROR3(1);
            ERROR_VEC(loc3+2) = ERROR3(2);
            ERROR_VEC(loc3+3) = ERROR3(3);

            if(max(abs(ERROR1)) > tol || max(abs(ERROR3)) > tol)
                concheck = 1;
            end
        end
    end   
end

end

