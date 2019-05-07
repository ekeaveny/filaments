classdef Filament < handle
% FILAMENT  Handle class for filaments.
%           N.B. This is a handle class and does NOT have a copy 
%           constructor. Assignment to a Filament type will cause both 
%           variables to be a reference to the same memory location.
    
properties

    N_w; % Number of segments comprising the filament.

    KB; % The elastic bending modulus of the filament.
    KT; % The elastic twisting modulus of the filament.

    DL; % The centreline distance between adjacent segments, Delta L.
    Length;  % Length of the filament, N_w*DL;
    
    StrainTwist; % The strain-twist vector encoding the intrinsic curvature and twist of the filament.

    weight_per_unit_length; % Weight per unit length.

    R; % The radii of the constituent particles.

    X; % Positions of the particles comprising the filament.
    Xm1; % time -1
    Xm2; % time -2

    X1; % Position of the first particle at the current and previous time-step.

    Q; % Orientation quaternions for the particles comprising the filament.

    U; % Lie algebra elements associated with the orientation of the constituent particles.

    V; % Velocities (both angular and translational) of the particles.

    F; % Forces and torques on the particles.

    Lambda; % The collection of Lagrange multipliers associated with the inextensibility of the filament.

    Lmat; % The lower-triangular part of the Jacobian L-U decomposition.
    Umat; % The upper-triangular part of the Jacobian L-U decomposition.

end

methods

    function obj = Filament(varargin)
    % FILAMENT(N_w, a, DL_factor)  constructs filament of N_w segments,
    %                              with thickness 2*a, and segment
    %                              separation a*DL_factor.
        
        if nargin==0
            % Empty constructor doesn't actually need to do anything,
            % it's only so that we can pre-allocate arrays to contain
            % Filaments.
        else
            if nargin==1
                a = 1;
                DL_factor = 2.2;
            elseif nargin==2
                a = varargin{2};
                DL_factor = 2.2;
            else
                a = varargin{2};
                DL_factor = varargin{3};
            end
                
            N = varargin{1};                
            obj.N_w     = N;               
            obj.R      = ones(1,N)*a;    
            obj.X      = zeros(3,N);
            obj.Xm1    = zeros(3,N);
            obj.Xm2    = zeros(3,N);
            obj.X1     = zeros(6,1); 
            obj.Q      = zeros(N,8);
            obj.U      = zeros(9,N);
            obj.V      = zeros(6,N);
            obj.F      = zeros(6,N);
            obj.DL     = DL_factor*obj.R(1);
            obj.Length = obj.DL*obj.N_w;
            obj.Lambda = zeros(3,N-1);
        end
    end

    function RobotArm(obj)
    % Filament.ROBOTARM()  Starting at X(1), construct a beam like a 
    %                      robot arm, by adding a segment a distance 
    %                      Filament.DL  away, with the angle 
    %                      dictated as the average of the two segments' 
    %                      tangents. The process is repeated like this.
    %                      You must call  Filament.InitialSetup  before
    %                      calling this.

        Xtemp = obj.X;
        Qtemp = obj.Q;
        dL = obj.DL;
        for i = 2:obj.N_w
            Xtemp(:,i) = Xtemp(:,i-1) + ...
                dL/2*(QuaternionRotation(Qtemp(i-1,1:4),[1;0;0]) + ...
                QuaternionRotation(Qtemp(i,1:4),[1;0;0]));
        end
        obj.X = Xtemp;   
    end

    function InitialSetup(obj,FirstSegPos,StrainTwist,B,Alpha,Beta,...
                                                 weight_per_unit_length,mu)
    % Filament.INITIALSETUP(FirstSegPos,StrainTwist,B,Alpha,Beta,...
    %                                            weight_per_unit_length,mu)  
    %     Set first segment at position FirstSegPos, then use RobotArm to 
    %     construct the rest of the segments. Properties of the filament 
    %     are as the names of the inputs suggest.
        
        obj.X(:,1) = FirstSegPos;
        obj.X1 = [FirstSegPos;FirstSegPos];

        % The orientation quaternions are initialised to the identity by 
        % default.
        obj.Q = [ones(obj.N_w,1),zeros(obj.N_w,3),ones(obj.N_w,1), ...
                                                         zeros(obj.N_w,3)];

        obj.RobotArm;
        obj.Xm1 = obj.X;
        obj.Xm2 = obj.X;
        obj.StrainTwist = StrainTwist;
        obj.weight_per_unit_length = weight_per_unit_length;
        obj.KB = weight_per_unit_length*(obj.Length)^3 / B;
        obj.KT = obj.KB;
    end


    function InitialGuess(obj)
    % Filament.INITIALGUESS()  rearranges linear interpolation as a guess 
    %                          for x^(j+1), i.e.
    %                          x^j = 0.5*( x^(j-1) + x^(j+1) ) .
    %                          Only does position of particle 1, because
    %                          RobotArm does the rest.
        
        % Guess position of particle 1.
        obj.X(:,1) = 2*obj.X1(1:3) - obj.X1(4:6); 
        Utemp = obj.U;
        Qtemp = obj.Q;

        for i=1:obj.N_w
            % Guess Lie algebra elements, and hence quaternions, of all
            % particles.
            Utemp(1:3,i) = 2*Utemp(4:6,i) - Utemp(7:9,i);
            Qtemp(i,1:4) = QuaternionProduct(qexp(Utemp(1:3,i)),...
                                                             Qtemp(i,5:8));
        end
        obj.U = Utemp;
        obj.Q = Qtemp;
    end

    function InternalForcesAndTorques(obj)
    % Filament.INTERNALFORCESANDTORQUES()  Place gravity, elastic and 
    %                                      constraint forces and torques 
    %                                      on the filament. They are stored
    %                                      in Filament.F

        N = obj.N_w;

        % Gravity (in z-direction)
        ForcesAndTorques = [zeros(2,N);...
                            -obj.weight_per_unit_length*obj.DL*ones(1,N);...
                            zeros(3,N)]; 

        % Elastic torques on segment 1
        ForcesAndTorques(4:6,1) = ForcesAndTorques(4:6,1) + ...
                                     obj.KB/(obj.N_w * obj.DL) * [1; 1; 1];

        Qtemp = obj.Q;
        Lam = obj.Lambda;
        dL = obj.DL;
        ST = obj.StrainTwist;
        Bend = obj.KB;
        Twist = obj.KT;

        for i=1:obj.N_w-1
            % Constraint forces and torques
            ForcesAndTorques(1:3,i) = ForcesAndTorques(1:3,i) - Lam(:,i);

            ForcesAndTorques(1:3,i+1) = ForcesAndTorques(1:3,i+1) ...
                                      + Lam(:,i);

            t = QuaternionRotation(Qtemp(i,1:4),[1;0;0]);

            ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) - ...
                0.5*dL*[t(2)*Lam(3,i) - t(3)*Lam(2,i);t(3)*Lam(1,i) ...
                                      - t(1)*Lam(3,i);t(1)*Lam(2,i) ...
                                      - t(2)*Lam(1,i)];

            t = QuaternionRotation(Qtemp(i+1,1:4),[1;0;0]);

            ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - ...
                0.5*dL*[t(2)*Lam(3,i) - t(3)*Lam(2,i);t(3)*Lam(1,i) ...
                      - t(1)*Lam(3,i);t(1)*Lam(2,i) - t(2)*Lam(1,i)];

            % Elastic torques
            q = MidpointQ(Qtemp(i,1:4),Qtemp(i+1,1:4));

            dqds = (Qtemp(i+1,1:4) - Qtemp(i,1:4))/dL;

            b = 2*QuaternionProduct([q(1),-q(2),-q(3),-q(4)],dqds);
            b = b(2:4);

            M = QuaternionRotation(q,[Twist * (b(1) - ST(3)); ...
                Bend * (b(2) - ST(1)); Bend * (b(3) - ST(2))]);

            ForcesAndTorques(4:6,i) = ForcesAndTorques(4:6,i) + M;

            ForcesAndTorques(4:6,i+1) = ForcesAndTorques(4:6,i+1) - M;
        end
        obj.F = ForcesAndTorques;
    end

    function com = CentreOfMass(obj)
    % Filament.CENTREOFMASS()  returns the centre of mass of the
    %                          filament.
        com = mean(obj.X,2);
    end

    function ConstructAndDecomposeJacobian(obj,dt,mu)
    % Filament.CONSTRUCTANDDECOMPOSEJACOBIAN(dt,mu)  constructs the
    %     approximate Jacobian for that filament. It places the
    %     decomposed parts into Filament.Lmat and Filament.Umat.
        [obj.Lmat,obj.Umat] = lu(approximate_jacobian(obj,dt,mu)); 
    end

    function out = InvertLocalBlock(obj,v)
    % Filament.INVERTLOCALBLOCK(v)  returns J_0\v, where J_0 is the
    %                               pre-computed, pre-LU-decomposed
    %                               approximate Jacobian for that 
    %                               filament.
        out = obj.Lmat\v;
        out = obj.Umat\out;
    end

    function ApplyUpdate(obj,u)
    % Filament.APPLYUPDATE(u)  Alg 2, Line 5. Does  X = X + u  for state
    %                          vector X.
        N = obj.N_w;
        obj.X(:,1) = obj.X(:,1) + u(1:3);
        Utemp = obj.U;
        Qtemp = obj.Q;
        Lam = obj.Lambda;
        for i=1:N
            Utemp(1:3,i) = Utemp(1:3,i) + u(3*i+1:3*(i+1));
            Qtemp(i,1:4) = QuaternionProduct(qexp(Utemp(1:3,i)),...
                                                             Qtemp(i,5:8));
            if i<N
                Lam(:,i) = Lam(:,i) + u(3*(N+i)+1:3*(N+i+1));
            end
        end
        obj.U = Utemp;
        obj.Q = Qtemp;
        obj.Lambda = Lam;
    end

    function EndOfStepUpdate(obj)
    % Filament.ENDOFSTEPUPDATE()  Step in time. Sets historical values of 
    %                             positions (etc.) to their new values.
        
        obj.X1(4:6) = obj.X1(1:3);
        obj.X1(1:3) = obj.X(:,1);

        obj.Xm2 = obj.Xm1;
        obj.Xm1 = obj.X;

        Utemp = obj.U;
        Qtemp = obj.Q;
        Utemp(7:9,:) = Utemp(4:6,:);
        Utemp(4:6,:) = Utemp(1:3,:);

        for i=1:obj.N_w
            Qtemp(i,5:8) = QuaternionProduct(qexp(Utemp(1:3,i)),...
                                                             Qtemp(i,5:8));
        end

        obj.U = Utemp;
        obj.Q = Qtemp;
    end
    
    function PrintToFile(obj,fid,t)
    % Filament.PRINTTOFILE(fid,t)  saves a string of filament segments 
    %                              positions, quaternions, velocities and 
    %                              forces to the file handle provided by  
    %                              fid .
    
        for j_w = 1:obj.N_w
            if j_w == obj.N_w
                L1 = 0;
                L2 = 0;
                L3 = 0;
            else
                L1 = obj.Lambda(1,j_w);
                L2 = obj.Lambda(1,j_w);
                L3 = obj.Lambda(1,j_w);                    
            end
            fprintf(fid, ['%.2f' repmat(' %.6f',1,22) '\n'], ...
            t, obj.X(1,j_w), obj.X(2,j_w), obj.X(3,j_w), ...
            obj.Q(j_w,1), obj.Q(j_w,2), obj.Q(j_w,3), obj.Q(j_w,4), ...
            obj.V(1,j_w), obj.V(2,j_w), obj.V(3,j_w), ...
            obj.V(4,j_w), obj.V(5,j_w), obj.V(6,j_w), ...
            obj.F(1,j_w), obj.F(2,j_w), obj.F(3,j_w), ...
            obj.F(4,j_w), obj.F(5,j_w), obj.F(6,j_w), ...
            L1, L2, L3);
        end
    end
        
end
    
end

