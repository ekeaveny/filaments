function main()
% MAIN  Supplementary code to 'Methods for suspensions of passive and
%       active filaments', https://arxiv.org/abs/1903.12609 ,
%       by SF Schoeller, AK Townsend, TA Westwood & EE Keaveny.
%       Visit https://github.com/ekeaveny/filaments/ to find contact 
%       details.
%       This version: 7 May 2019
%
%   This code demonstrates the use of the method in simulating a
%   single flexible filament falling under gravity.
%
%   It uses the 'EJBb' version of Broyden's method (Algorithm 2 in the
%   paper) with a reduced 'robot arm' system of nonlinear equations.
%
%   To use, just run this script.


% Setup
save_to_file = true;
graphics = true;
plot_step = 5;                % Plot every n timesteps
save_step = 5;                % Save data to file every n timesteps

% Filament data
a = 1;                        % segment 'radius' (half filament width)
N_sw = 1;                     % number of filaments
N_w = 30;                     % number of segments in filament
Np = N_sw*N_w;                % total number of segments
max_broyden_steps = 3*Np;

B = 1e3;                      % dimensionless elasto-gravitational number B
weight_per_unit_length = 1e0; % weight per unit length W

steps_per_unit_time = 300;    % for sedimenting, unit time T = L^2 mu / F
num_settling_times = 2;       % number of settling (unit) times
concheck_tol = 1e-4;          % Broyden's tolerance

DL_factor = 2.2;              % distance between segment centres, Delta L,
                              % is DL = DL_factor*a.
mu = 1;                       % fluid viscosity

filename = ['out-'  datestr(now,'yyyymmdd-HHMMSS') '-Nsw' num2str(N_sw) ...
            '-Nw' num2str(N_w) '-B' num2str(B) '-RPY']; % used for data
                                                        % file & video file
% Set up segment positions.
% (Note it does it from N_sw to 1 so that the number of filaments is only
% changed once.)
for i = N_sw:-1:1
    % Initialise filament object
    Filaments(i) = Filament(N_w,a,DL_factor);
    % Position filament
    Filaments(i).InitialSetup([0;10*(i-1);0],[0,0,0],B,0,0, ... 
                                                weight_per_unit_length,mu);
    N(i) = 6*Filaments(i).N_w;
end

N = [0,cumsum(N)];
Nbroy = N(end);

% Time
unit_time = Filaments(1).Length*mu/weight_per_unit_length; % 1 settling 
                                                           % time, T
TOTAL_STEPS = num_settling_times*steps_per_unit_time;
dt = unit_time/steps_per_unit_time;
t = 0;
plot_now = plot_step - 1;
save_now = save_step - 1;

% Time and iteration counts
frame_time = zeros(TOTAL_STEPS,1);
iters = zeros(TOTAL_STEPS,1);       % Number of Broyden's iterations
running_total_count = 0;            % For average number of Broyden's iters

for nt = 1:TOTAL_STEPS
    iter = 0;
    
    p_broy = max_broyden_steps + 1;
    Cmat = zeros(Nbroy,p_broy); % c and d vectors from Alg 2, Line 7. Their 
    Dmat = zeros(Nbroy,p_broy); % value at each iteration is stored.
    
    % Stop if broken
    if isnan(Filaments(1).X(1))
        keyboard
        continue
    end    
    
    % Screen output
    fprintf('\n')
    if mod(nt,20) == 0
        fprintf(['[' filename ': rEJBb, B=' num2str(B) ', RPY, Nsw=' ...
                 num2str(N_sw) ', Nw=' num2str(N_w) ']\n' ])
    end
    length_of_TOTAL_STEPS = max(ceil(log10(abs(TOTAL_STEPS))),1);
    fprintf([ 'B=' num2str(B) '   ' ...
              'timestep: ' ...
              sprintf(['%' num2str(length_of_TOTAL_STEPS) '.f'],nt) ...
              '/' num2str(TOTAL_STEPS) ' ' ])
    frame_start = tic;    
    
    % Aim of this is to update positions
    for i=1:N_sw
        % Guess new position of first segment in filament, and
        % new quaternions for all segments in filament
        Filaments(i).InitialGuess;
        % use RobotArm to guess rest of positions
        Filaments(i).RobotArm;
    end
    
    % Find f(X_k) and place into ERROR_VECk.
    % If ||ERROR_VECk|| < concheck_tol (= epsilon in Alg 2, Line 4),
    % then concheck = 0. Else, 1.    
    [concheck,ERROR_VECk] = F(Filaments,concheck_tol); 
    
    % Find approximate Jacobian J_0. Since this is found in filament-sized
    % blocks, it is stored as property of the filament objects.
    % For convenience, it's also LU-decomposed at this stage.
    for i=1:N_sw
        Filaments(i).ConstructAndDecomposeJacobian(dt,mu);
    end
    
    % Find J_0^{-1} f(X_k)  (from Alg 2, Line 5)
    J0invERROR_VECk = blockwise_backslash_jacobian(Filaments, ERROR_VECk);    
     
    num_broydens_steps_required = 0;
    while (concheck == 1) % Alg 2, Line 4
        % Alg 2, Line 5. DeltaX is Delta X in paper.
        DeltaX = -apply_inverse_jacobian(J0invERROR_VECk, Cmat, Dmat,...
                                         ERROR_VECk, iter);

        % Update the positions and lambdas
        for i=1:N_sw
            Filaments(i).ApplyUpdate(DeltaX(N(i)+1 : N(i+1)));
            Filaments(i).RobotArm;
        end
        
        % Check to see if the new state is an acceptable solution:
        % ERROR_VECk1 = f(X_(k+1))        
        [concheck,ERROR_VECk1] = F(Filaments,concheck_tol);     
               
        iter = iter + 1;
        
        % (remaining lines are Alg 2, Line 7)
        y_vec = ERROR_VECk1 - ERROR_VECk;
        
        J0invERROR_VECk1 = blockwise_backslash_jacobian(Filaments, ...
                                                              ERROR_VECk1);
        
        y_vec_sq = y_vec'*y_vec;
        Cmat(:,iter) = -apply_inverse_jacobian(J0invERROR_VECk1, ...
                                            Cmat, Dmat, ERROR_VECk1, iter);
        Dmat(:,iter) = y_vec/y_vec_sq;
        ERROR_VECk = ERROR_VECk1;
        J0invERROR_VECk = J0invERROR_VECk1;
        
        % Shout if the iteration count has got a bit high
        if iter == 100
            keyboard
            continue
        end

        % If the number of iterations maxes out, proceed to next timestep
        % anyway and see what happens (but flag it with a *)
        if (iter > max_broyden_steps)
            fprintf(' *');
            concheck = 0;
        end

        num_broydens_steps_required = num_broydens_steps_required + 1;
        running_total_count = running_total_count + 1;        
    end
    
    % Step in time, step in time. Never need a reason, never need a rhyme..
    t = t + dt;    
    for i=1:N_sw
        Filaments(i).EndOfStepUpdate;
    end
    
    % Plot and save
    plot_now = plot_now + 1;
    save_now = save_now + 1; 
    if(save_now == save_step && save_to_file)
        fid = fopen(['output/' filename '.dat'], 'a');
        if nt == 1
            fprintf(fid, 'dt, B, Nsw (RPY)\n');
            fprintf(fid, '%.6e %.6e %.f\n\n', dt, B, N_sw);
        end
        for j_sw = 1:N_sw
            Filaments(j_sw).PrintToFile(fid,t);
        end
        fprintf(fid,'\n');
        fclose(fid);
        clf;
    end    
    
    if (plot_now == plot_step && graphics)
        com = Filaments(1).CentreOfMass;
        L = Filaments(1).Length;
        for i_sw = 1:N_sw
            x = Filaments(i_sw).X;
            plot3(x(1,:)/L,x(2,:)/L,x(3,:)/L,'-','LineWidth',10);
            hold on;
        end
        
        % Title
        title(['nt='  num2str(nt)  ', dt='  num2str(dt) ...
               ', B='  num2str(B)  ', N_{sw}='  num2str(N_sw) ...
               ])        
        
        hold off
        pbaspect([1 1 1])
        xlim([com(1)/L-0.5,com(1)/L+0.5]);
        ylim([com(2)/L-0.5,com(2)/L+0.5]); 
        zlim([com(3)/L-0.5,com(3)/L+0.5]); 
        xlabel('(x-x_{COM})/L');
        ylabel('(y-y_{COM})/L');
        zlabel('(z-z_{COM})/L');

        pause(0.01);
    end 
    
    if plot_now == plot_step
        plot_now = 0;
    end
    if save_now == save_step
        save_now = 0;
    end
    
    frame_time(nt) = toc(frame_start);
    iters(nt) = iter;

    fprintf(['[' format_time(frame_time(nt)) '|' ...
            format_time(mean(frame_time(1:nt))*(TOTAL_STEPS-nt)) ...
            '-][#Broy steps: '  num2str(num_broydens_steps_required) ...
            '|Avg: '  num2str(round(running_total_count/nt))  ']'])    
   
end

disp('')
disp('Run finished')
disp(['Total time:' format_time(sum(frame_time))])


function [concheck_local,ERROR_VECk1_local] = F(Filaments,tol)
% F  places forces and torques on the segments, calculates the resultant
%    velocities and angular velocities, and forms the error vector f(X*).
%    Then checks convergence. For details, see docstrings of functions
%    within.

    for ii=1:N_sw
        Filaments(ii).InternalForcesAndTorques();
    end
    collision_barrier(Filaments);
    
    RPY(Filaments,mu);
    
    % Check convergence between x_(n+1) and x_n, and also check the
    % constraint. concheck = 0 if all fine, 1 otherwise. The error vectors
    % are all compiled into ERROR_VECk1_local.    
    [concheck_local,ERROR_VECk1_local] = constraint_check_robot_arm(...
                                                      Filaments,dt,nt,tol);
end

end

