function [KKT, rhs] = createKKTsparse(H,g,A,b)
% KKTsparse   Create a sparse KKT system
% 
% Creates a KKT sparse matrix (using the built in sparse matrix in Matlab)
% and matching right hans side from a EQP formulation. Optionally makes the KKT matrix sparse.
%
% Inputs:
%   H   : a positive semi-definite Hessian matrix of dimension (n,n) defining the quadratic terms of the objective function
%   g   : a (n,) dimensional vector defining the linear part of the objective function
%   A   : a (n,m) dimensional matrix giving the rhs of the equality constrains
%   b   : a (m,) dimensional vector defining the rhs of the equality constrains
%
% Outputs:
%   KKT : a KKT sparse matrix with dimensions (2n, n+m)
%   rhs : a vector of dimension (n+m,) the right hand side of the KKT system
%
%%

[n,m] = size(A);

[KKT, rhs] = createKKT(H,g,A,b);
KKT = sparse(KKT);

end