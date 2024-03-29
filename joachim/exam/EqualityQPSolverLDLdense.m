function [x, lambda, time] = EqualityQPSolverLDLdense(H,g,A,b)
% EqualityQPSolverLDLdense  EQP solver using LDL-factorization, uses a dense KKT system.
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

[n,m] = size(A);

% Create a the KKT system
[KKT, rhs] = createKKT(H,g,A,b);

% Fatorize KKT matrix
tstart = cputime;
[L,D,p] = ldl(KKT,'lower','vector');
time = cputime-tstart;

% Solve for x and lagrange multipliers
z = zeros(n+m, 1);
z(p) = L' \ (D \ (L \ rhs(p)));
x = z(1:n);
lambda = z(n+1:n+m);


