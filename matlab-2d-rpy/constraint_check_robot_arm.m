function [concheck, ERROR_VEC] = constraint_check_robot_arm(...
                                  X_S, Y_S, THETA_S, X, Y, THETA, ...
                                  X_T, Y_T, THETA_T, VX, VY, OMEGZ, ...
                                  DL, dt, nt, SW_IND, tol)
% CONSTRAINT_CHECK_ROBOT_ARM  Creates the error vector for the Broyden's
%                             step and finds whether it is less than a 
%                             specified tolerance.
%
%   This corresponds to the system f(X*) at equation (41) in the paper.
%
%   [concheck, ERROR_VEC] = constraint_check_robot_arm(...
%                                  X_S, Y_S, THETA_S, X, Y, THETA, ...
%                                  X_T, Y_T, THETA_T, VX, VY, OMEGZ, ...
%                                  DL, dt, nt, SW_IND, tol)
%   takes in calculated next timestep positions/orientations (X_S etc.),
%                    current timestep positions/orientations (X   etc.),
%               and previous timestep positions/orientations (X_T etc.),
%               the calculated velocities/angular velocities (VX  etc.),
%   segment separation DL, timestep dt, timestep number nt, 
%   filament indexing matrix SW_IND, and tolerance (epsilon in paper) tol.
%                                                                  
%   It returns ERROR_VEC (f(X) in the paper) and whether ||f(X)|| < tol.
%   concheck = 0 means ||f(X)|| < tol. Otherwise, concheck = 1.

Np = size(X, 1);

TX_S = cos(THETA_S);
TY_S = sin(THETA_S);

N_sw = size(SW_IND,1);
N_w = size(SW_IND,2);

concheck = 0;

ERROR_VEC = zeros(3*Np,1);

c_1 = 4.0/3.0;
c_2 = 1.0/3.0;
c_4 = 2.0/3.0;

if(nt == 1) % Euler
    for j_sw = 1:N_sw
        first_segment_index = SW_IND(j_sw,1);
        TX_S_sum = 0;
        TY_S_sum = 0;
        for j_w = 1:N_w
            j = SW_IND(j_sw,j_w); % segment number j
            TX_S_sum = TX_S_sum + TX_S(j); % } sum[TX_S(1:j)] but general-
            TY_S_sum = TY_S_sum + TY_S(j); % } ised for multiple filaments
            ERROR1 = X_S(first_segment_index) + DL*TX_S_sum ...
                   - DL/2*(TX_S(first_segment_index) + TX_S(j)) ...
                   - X(j) - dt*VX(j);
            ERROR2 = Y_S(first_segment_index) + DL*TY_S_sum ...
                   - DL/2*(TY_S(first_segment_index) + TY_S(j)) ...
                   - Y(j) - dt*VY(j);           
            
            ERROR_VEC(j) = ERROR1;
            ERROR_VEC(Np + j) = ERROR2;
            if(abs(ERROR1) > tol || abs(ERROR2) > tol)
                concheck = 1;
            end        
        end
    end
        
    for j = 1:Np
        ERROR3 = THETA_S(j) - (THETA(j) + dt*OMEGZ(j));
        ERROR_VEC(2*Np + j) = ERROR3;
        if(abs(ERROR3) > tol)
            concheck = 1;
        end         
    end
                  
else % BDF2
    for j_sw = 1:N_sw
        first_segment_index = SW_IND(j_sw,1);
        TX_S_sum = 0;
        TY_S_sum = 0;        
        for j_w = 1:N_w
            j = SW_IND(j_sw,j_w); % segment number j
            TX_S_sum = TX_S_sum + TX_S(j); % } sum[TX_S(1:j)] but general-
            TY_S_sum = TY_S_sum + TY_S(j); % } ised for multiple filaments          
            ERROR1 = X_S(first_segment_index) + DL*TX_S_sum ...
                   - DL/2*(TX_S(first_segment_index) + TX_S(j)) ...
                   - c_1*X(j) + c_2*X_T(j) - dt* c_4*VX(j);
            ERROR2 = Y_S(first_segment_index) + DL*TY_S_sum ...
                   - DL/2*(TY_S(first_segment_index) + TY_S(j)) ...
                   - c_1*Y(j) + c_2*Y_T(j) - dt* c_4*VY(j);
            ERROR_VEC(j) = ERROR1;
            ERROR_VEC(Np + j) = ERROR2;   
            if(abs(ERROR1) > tol || abs(ERROR2) > tol)
                concheck = 1;
            end        
        end
   end
   for j=1:Np
        ERROR3 = THETA_S(j) ...
               - (c_1*THETA(j) - c_2*THETA_T(j) + dt* c_4*OMEGZ(j));
        ERROR_VEC(2*Np + j) = ERROR3;
        if(abs(ERROR3) > tol)
            concheck = 1;
        end            
   end
   
   
end

