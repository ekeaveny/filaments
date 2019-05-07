function [A] = rcross(x)
% RCROSS(x)  Given a column vector x, returns the matrix A such that 
%            A*v = cross(v,x);

A = [0,x(3),-x(2);-x(3),0,x(1);x(2),-x(1),0];

end

