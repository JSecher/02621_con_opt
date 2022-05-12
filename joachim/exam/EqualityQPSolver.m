function [x, lambda, time] = EqualityQPSolver(H, g, A, b, solver)
% EqualityQPSolver  Interface for variuos different EQP solvers.
%
% Solves EQP problems on the form:
%
%   \min_{x} f(x) = x^{T}.H.x + g^{T}.x
%    s.t.    A.x  = b  
%
% With options for using solvers with different facotrizations. 
% 
% Inputs:
%   H       : a positive semi-definite Hessian matrix of dimension (n,n) defining the quadratic terms of the objective function
%   g       : a (n,) dimensional vector defining the linear part of the objective function
%   A       : a (n,m) dimensional matrix giving the rhs of the equality constrains
%   b       : a (m,) dimensional vector defining the rhs of the equality constrains
%   solver  : a string giving the solver to be used, availible obtions are: 'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace'
%
% Outputs:
%   x       : a (n,) dimensional vector of found solution
%   lambda  : a (m,) dimensional vector of lagrange multipliers
%   time    : a float giing the CPU time in seconds on the factorization
%
% See also EqualityQPSolverLDLdense, EqualityQPSolverLDLsparse, EqualityQPSolverLUdense, EqualityQPSolverLUsparse, EqualityQPSolverRangeSpace, EqualityQPSolverNullSpace
%%

solver_match = lower(solver);
switch solver_match
    case "ludense"
        [x, lambda,time] = EqualityQPSolverLUdense(H,g,A,b);
    case "lusparse"
        [x, lambda,time] = EqualityQPSolverLUsparse(H,g,A,b);
    case "ldldense"
        [x, lambda,time] = EqualityQPSolverLDLdense(H,g,A,b);
    case "ldlsparse"
        [x, lambda,time] = EqualityQPSolverLDLsparse(H,g,A,b);
    case "rangespace"
        [x, lambda,time] = EqualityQPSolverRangeSpace(H,g,A,b);
    case "nullspace"
        [x, lambda,time] = EqualityQPSolverNullSpace(H,g,A,b);
    otherwise
        error("Unknown solver choice: '" + solver + "'. Availible options are: 'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace' and 'NullSpace'.")
end

