function [X_S,Y_S] = robot_arm(X_S,Y_S,THETA_S,SW_IND,DL)
% ROBOT_ARM  Starting at [X_S(1),Y_S(1)], construct a beam like a robot 
%            arm, by adding a segment a distance DL away, with the angle 
%            dictated as the average of the two segments' tangents. The 
%            process is repeated like this.
%
%   [X_S,Y_S] = robot_arm(X_S,Y_S,THETA_S,SW_IND,DL)
%
%   The tangent vector [TX_S,TY_S] is generated from THETA_S. Which segment
%   is in which filament is conveyed through SW_IND.
% 
%   You should specify [X_S(1),Y_S(1)] for each filament before calling
%   this function.

N_sw = size(SW_IND,1);
    N_w = size(SW_IND,2);
    TX_S_loc = cos(THETA_S);
    TY_S_loc = sin(THETA_S);
    for k_sw = 1:N_sw
        for k_w = 2:N_w
            iii = SW_IND(k_sw,k_w);
            iii_prev = SW_IND(k_sw,k_w-1);
            X_S(iii) = X_S(iii_prev) + DL/2*(TX_S_loc(iii_prev) + TX_S_loc(iii));
            Y_S(iii) = Y_S(iii_prev) + DL/2*(TY_S_loc(iii_prev) + TY_S_loc(iii));
        end
    end    
end