function [x, lambda, time] = EqualityQPSolverNullSpace(H,g,A,b)
% EqualityQPSolverNullSpace EQP solver using the  Null Space method
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

% Factorize A matrix using QR decomposition
tstart = cputime; 
[Q,Rbar] = qr(A, 'vector');
time = cputime-tstart;

if ~ m == size(Rbar, 2)
    error("A does not have full column rank")
end

% Solve for x and lagrange multipliers    
Y = Q(:,1:m);      % Q Range
Z = Q(:,m+1:n);    % Q Null
R = Rbar(1:m,1:m);
x_Y = R' \ b;       % Solve R' x_Y = b

% Solve: (Z' H Z)x_Z = −Q′2(HQ1 x_Y + g) i.e. solve for x_Z
ZT = Z';           % Precompute tranpose Z for use multiple times
L = chol(ZT*H*Z);  % Chol factorize the reduced-Hessian matrix

x_Z = L' \ (-ZT * (H*(Y*x_Y)+g));  % Solve for x_Z
x_Z = L \ x_Z;                       

x = Y*x_Y+Z*x_Z;              % Compute x = Y x_Y + Z x_Z
lambda = R \ Y'*(H*x + g);     % Solve: R lambda = Q1′(Hx + g)

end
