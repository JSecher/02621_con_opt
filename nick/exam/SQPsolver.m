function [x, obj, lambda, output] = SQPsolver(objfun,confun,xlower,xupper,clower,cupper,x0,solver)
% SQPSolver  Interface for variuos SQP solvers.
%
% Solves NLP's on the form  
%
%   \min_{x} f(x) 
%    s.t.    cl <= c(x) <= cu  
%            xl <= x <= xu
%
% With options for using solvers with different facotrizations. 
% 
% Inputs:
%   objfun  : Objective function, should take a vector x of size size(n),
%           : and return [function value at x, gradients at x]
%   confun  : Constraint function, should take a vector x of size size(n)
%           : and reutrn [constraint function value at x, constraint gradients at x]
%   xlower  : (n,)-dim vector of lower bounds for x
%   xupper  : (n,)-dim vector of upper bounds for x
%   clower  : (m,)-dim vector of lower bounds for x
%   cupper  : (m,)-dim vector of lower bounds for x
%   x0      : (n,)-dim vector of inital point for x
%   solver  : string with the solver of choice, options are "bfgs", "line",
%             and "trust"
%
% Outputs:
%   x       : a (n,) dimensional vector of found solution
%   obj     : a float giving the objective value
%   lambda  : a (m,) dimensional vector of lagrange multipliers
%   output  : a object with performance and iteration information
%
% See also SQPsolverBFGS
%%

solver_match = lower(solver);
switch solver_match
    case "bfgs"
        [x, lambda, output] = SQPsolverBFGS(objfun,confun,xlower,xupper,clower,cupper,x0);
    case "line"
        error("Line search is not implemented")
    case "trust"
        error("Trust region is not implemented")
    otherwise
        error("Unknown solver choice: '" + solver + "'. Availible options are: 'bfgs', 'line', and 'trust'.")
end
obj = objfun(x);

