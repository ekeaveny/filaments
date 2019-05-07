function [out] = MidpointQ(q1,q2)
% MIDPOINTQ   Given 2 unit quaternions, q1 and q2, this function outputs 
%             the quaternion 'half way' between them in a rotational sense.

% First, we calculate the rotation mapping q1 to q2.

% Conjugate of q1.
q1_conj = [q1(1),-q1(2),-q1(3),-q1(4)]; 

% q_rot is the quaternion mapping q1 to q2.
q_rot = QuaternionProduct(q2,q1_conj); 

% Produce the half rotation.
half_rot = unit_q_sqrt(q_rot); 

% Produce the midpoint.
out = QuaternionProduct(half_rot,q1); 

end

