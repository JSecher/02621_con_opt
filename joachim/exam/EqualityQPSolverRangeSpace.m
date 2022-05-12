function [x, lambda, time] = EqualityQPSolverRangeSpace(H,g,A,b)
% EqualityQPSolverNullSpace EQP solver using the Range Space method
%
% Inputs:
%   H       : a positive semi-definite Hessian matrix of dimension (n,n) defining the quadratic terms of the objective function
%   g       : a (n,) dimensional vector defining the linear part of the objective function
%   A       : a (n,m) dimensional matrix giving the rhs of the equality constrains
%   b       : a (m,) dimensional vector defining the rhs of the equality constrains
%
% Outputs:
%   x       : a (n,) dimensional vector of found solution
%   lambda  : a (m,) dimensional vector of lagrange multipliers
%   time    : a float giing the CPU time in seconds on the factorization
%
% See also EqualityQPSolver
%%

% Factorize H matrix
tstart = cputime;
L = chol(H); 
time = cputime-tstart;

% Solve for x and lagrange multipliers    
LT = L';  % Precompute to avoid doing it twice
HA = L \ (LT \ A);
Hg = L \ (LT \ g);
lambda = (A'*HA) \ (b+A'*Hg);
x = (HA*lambda)-Hg;

end