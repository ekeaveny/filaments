function [J0] = approximate_jacobian(THETA, LAMBDA1, LAMBDA2, ...
                  drag_coeff, vistor_coeff, dt, DL, KB_div_DL_array, ...
                  SW_IND)
% APPROXIMATE_JACOBIAN  Creates a condensed blockwise Jacobian based on
%                       analytic expressions given in Appendix B of the
%                       paper.
%
%   The Jacobian is of the form [A1 A2 A3 ...], which is a compact
%   version of the form J_0 = [A1  O  O ..
%                               O A2  O ..
%                               O  O A3 ..
%                               :  :  :   ].
%
%   Each block A1, A2, etc., is of size (3*N_w, 3*N_w), and is the
%   Jacobian of each individual filament. Non-block diagonal elements
%   of the Jacobian are taken care of outside this function by Broyden's
%   method.
%
%   [J0] = approximate_jacobian(THETA, LAMBDA1, LAMBDA2, ...
%                               drag_coeff, vistor_coeff, dt, DL, ...
%                               KB_div_DL_array, SW_IND)
%   takes in angles THETA, constraints [LAMBDA1,LAMBDA2], drag and viscous
%   torque coefficients, timestep dt, segment separation vector DL, segment
%   indexing matrix SW_IND, and KB_div_DL_array, which allows for
%   different bending moduli for each filament.
%
%   J0 is the approximate Jacobian, J_0.

Np = size(THETA, 1);

N_sw = size(SW_IND,1);
N_w = size(SW_IND,2);
N_lam_w = N_w - 1;

J0 = zeros(3*N_w, 3*Np);

c1 = dt/drag_coeff;
c2 = dt/vistor_coeff;

for i_sw = 1:N_sw

    num_filament_types = size(KB_div_DL_array,2);
    k_index = floor((i_sw-1)/N_sw*num_filament_types) + 1;
    KB_div_DL = KB_div_DL_array(k_index);
    
    % Note order of the Jacobian produced is particle-sized
    %                                        [11 13 12
    %                                         31 33 32
    %                                         21 23 22]

    % [df1/dY1 = J11 and df3/dY1 = J31]
    for i_w = 1:N_w
        J = 3*(i_w-1)+1;                  % = j mod N_w
        k = 3*((i_sw-1)*N_w+  1-1)+1;
        J0(J,   k) = 1;                % 1s for Xs
        J0(J+1,k+1) = 1;               % 1s for Ys
    end

    % [df1/dlambda = J13 and df3/dlambda = J33]
    for i_w = 1:N_w-1
        j = 3*((i_sw-1)*N_w+i_w-1)+1;
        J = 3*(i_w-1)+1;                  % = j mod N_w
        J0(J,   j+3) = c1;             % c1s for Xs
        J0(J+3, j+3) = -c1;            % c1s for Xs
        J0(J+1, j+4) = c1;             % c1s for Ys
        J0(J+4, j+4) = -c1;            % c1s for Ys
    end

    % [df2/dtheta = J22]
    for i_w = 1:N_w
        ip = SW_IND(i_sw,i_w);
        j = 3*((i_sw-1)*N_w+i_w-1)+1;
        J = 3*(i_w-1)+1;                  % = j mod N_w
        if(i_w == 1)
            jp = SW_IND(i_sw,i_w + 1);
            ilam = (i_sw - 1)*N_lam_w + i_w;

            J0(J+ 2,j+2) = 1.0;
            J0(J+2,j+2) = J0(J+2,j+2) ...
                           + c2*KB_div_DL*(sin(THETA(ip))*sin(THETA(jp))...
                           + cos(THETA(ip))*cos(THETA(jp)));
            J0(J+2,j+2) = J0(J+2,j+2) ...
                           - c2*0.5*DL*(LAMBDA1(ilam)*cos(THETA(ip))...
                           + LAMBDA2(ilam)*sin(THETA(ip)));
            J0(J+2,j+5) = -c2*KB_div_DL*(sin(THETA(ip))*sin(THETA(jp))...
                           + cos(THETA(ip))*cos(THETA(jp)));
        elseif(i_w == N_w)
            mp = SW_IND(i_sw,i_w - 1);
            ilamm1 = (i_sw - 1)*N_lam_w + i_w - 1;

            J0(J+2,j+2) = 1.0;
            J0(J+2,j+2) = J0(J+2,j+2) ...
                           + c2*KB_div_DL*(sin(THETA(ip))*sin(THETA(mp))...
                           + cos(THETA(ip))*cos(THETA(mp)));
            J0(J+2,j+2) = J0(J+2,j+2) ...
                           - c2*0.5*DL*(LAMBDA1(ilamm1)*cos(THETA(ip))...
                           + LAMBDA2(ilamm1)*sin(THETA(ip)));
            J0(J+2,j-1) = -c2*KB_div_DL*(sin(THETA(ip))*sin(THETA(mp))...
                           + cos(THETA(ip))*cos(THETA(mp)));
        else
            jp = SW_IND(i_sw,i_w + 1);
            mp = SW_IND(i_sw,i_w - 1);
            ilam = (i_sw - 1)*N_lam_w + i_w;
            ilamm1 = (i_sw - 1)*N_lam_w + i_w - 1;

            J0(J+2,j+2) = 1.0;
            J0(J+2,j+2) = J0(J+2,j+2) ...
                           + c2*KB_div_DL*(sin(THETA(ip))*sin(THETA(jp)) ...
                           + cos(THETA(ip))*cos(THETA(jp)));
            J0(J+2,j+2) = J0(J+2,j+2) ...
                           + c2*KB_div_DL*(sin(THETA(ip))*sin(THETA(mp)) ...
                           + cos(THETA(ip))*cos(THETA(mp)));
            J0(J+2,j+2) = J0(J+2,j+2) ...
                           - c2*0.5*DL*(LAMBDA1(ilam)*cos(THETA(ip)) ...
                           + LAMBDA2(ilam)*sin(THETA(ip)));
            J0(J+2,j+2) = J0(J+2,j+2) ...
                           - c2*0.5*DL*(LAMBDA1(ilamm1)*cos(THETA(ip)) ...
                           + LAMBDA2(ilamm1)*sin(THETA(ip)));

            J0(J+2,j+5) = -c2*KB_div_DL*(sin(THETA(ip))*sin(THETA(jp))...
                           + cos(THETA(ip))*cos(THETA(jp)));
            J0(J+2,j-1) = -c2*KB_div_DL*(sin(THETA(ip))*sin(THETA(mp))...
                           + cos(THETA(ip))*cos(THETA(mp)));
        end
    end

    % [df2/dLambda = J23]
    for i_w = 1:N_w-1
        ip = SW_IND(i_sw,i_w);
        jp = SW_IND(i_sw,i_w + 1);
        j = 3*((i_sw-1)*N_w+i_w-1)+1;
        J = 3*(i_w-1)+1;                  % = j mod N_w

        J0(J+2,j+3) = -c2*0.5*DL*sin(THETA(ip));
        J0(J+5,j+3) = -c2*0.5*DL*sin(THETA(jp));

        J0(J+2,j+4) = c2*0.5*DL*cos(THETA(ip));
        J0(J+5,j+4) = c2*0.5*DL*cos(THETA(jp));
    end

    % [df3/dtheta = J32]
    for i_w = 1:N_w-1
        i1 = SW_IND(i_sw,1);
        i2top = SW_IND(i_sw,2:i_w);
        ip1 = SW_IND(i_sw,i_w+1);

        j = 3*((i_sw-1)*N_w+i_w-1)+1;
        J = 3*(i_w-1)+1;                  % = j mod N_w
        k = 3*((i_sw-1)*N_w+  1-1)+1;

        % X
        J0(J+3,k+2) = -0.5*DL*sin(THETA(i1));
        J0(J+3,j+5) = -0.5*DL*sin(THETA(ip1));
        J0(J+3,k+5:3:j+2) = -2*0.5*DL*sin(THETA(i2top));
        % Y
        J0(J+4,k+2) = 0.5*DL*cos(THETA(i1));
        J0(J+4,j+5) = 0.5*DL*cos(THETA(ip1));
        J0(J+4,k+5:3:j+2) = 2*0.5*DL*cos(THETA(i2top));
    end
end

end
