function [KKT, rhs] = createKKT(H,g,A,b)
% createKKT   Create a KKT system
%
% Creates a KKT matrix and matching right hans side
% from a EQP formulation. Optionally makes the KKT matrix sparse.
%
% Inputs:
%   H   : a positive semi-definite Hessian matrix of dimension (n,n) defining the quadratic terms of the objective function
%   g   : a (n,) dimensional vector defining the linear part of the objective function
%   A   : a (n,m) dimensional matrix giving the rhs of the equality constrains
%   b   : a (m,) dimensional vector defining the rhs of the equality constrains
%
% Outputs:
%   KKT : a KKT matrix with dimensions (n+m, n+m)
%   rhs : a vector of dimension (n+m,) the right hand side of the KKT system
%
%%

[n,m] = size(A);

KKT = [H, -A; -A', zeros(m)];
rhs = -[g; b];
end