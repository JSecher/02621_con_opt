function [H,g,A,b] = recycleEQP(n, uhat, d0)
% randEQP   Generate a EQP instance of recycle problem
%
% Inputs:
%   n       : a positive integer size of problem, must be greater or equal
%             to 3
%   uhat    : a float, parameter for the problem
%   scale   : a float, parameter for the problem
%
% Outputs:
%   H       : a symmetric positive semi-definite Hessian matrix of dimension (n,n) defining the quadratic terms of the objective function
%   g       : a (n,) dimensional vector defining the linear part of the objective function
%   A       : a (n,m) dimensional matrix giving the rhs of the equality constrains
%   b       : a (m,) dimensional vector defining the rhs of the equality constrains
%
%%

if n < 3
    error("n must be larger than 3")
end

H = eye(n+1);           % Generate Hessian
g = -uhat*ones(n+1,1);  % Set target
% Generate full A
A = -1*eye(n+1,n);      
A(1:(end-1), 2:end) = A(1:(end-1), 2:end) + eye(n,n-1);
A(end, end) = -1;
A(end-1, 1) = 1;
b = [-d0; zeros(n-1,1)];  % Generate rhs for constrainsts
    
