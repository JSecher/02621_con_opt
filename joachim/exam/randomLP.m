function [g,A,b,x,lam] = randomLP(n,m,scale)
% randEQP   Generate a random convex EQP
%
% Inputs:
%   n       : a positive integer giving number of variabels
%   m       : a positive integer giving number of constraints
%   scale   : a float which simply scales the values of the problem
%             (optional). Default value: 1
%
% Outputs:
%   H       : a symmetric positive semi-definite Hessian matrix of dimension (n,n) defining the quadratic terms of the objective function
%   g       : a (n,) dimensional vector defining the linear part of the objective function
%   A       : a (n,m) dimensional matrix giving the rhs of the equality constrains
%   b       : a (m,) dimensional vector defining the rhs of the equality constrains
%   x       : a (n,) dimensional vector giving the solution
%   lam     : a (m,) dimensional vector giving the lagrange multiplier
%
%%

if nargin < 3
    scale = 1;
end

H = n*eye(n);

% Generate A matrix with full rank
% Check if it is, retry if not
% (While loops and check is probably not nessesary, just 'in case')
maxrank = m;
rankA = 0;
triesA = 0;
while rankA < maxrank
    triesA = triesA + 1;
    if (triesA > 100)
        error("RuntimeError: Failed to generate valid A matrix in 100 tries")
    end
    A = scale*rand(n,m);
    rankA = rank(A);
end   

% Generate a random solution to the system
x = rand(n,1);
lam = rand(m,1);

% Generate matching g and b vectors from using KKT system
KKT = [H -A; -A' zeros(m)];
sol = [x; lam];
rhs = -(KKT*sol);
g = rhs(1:n);
b = rhs(n+1:end);

end







