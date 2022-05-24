function [x, lambda, iter] = SQPsolverBFGS(objfun,confun,xlower,xupper,clower,cupper,x0)
% SQPsolverBFGS SQP solvers using dampend BFGS approximation.
%
% Solves NLP's on the form  
%
%   \min_{x} f(x) 
%    s.t.    0 <= c(x) <= b  
%
% With options for using solvers with different facotrizations. 
% 
% Inputs:
%   objfun  : Objective function, should take a vector x of size size(n)
%   confun  : Objective function, should take a vector x of size size(n)
%   xlower  : (n,)-dim vector of lower bounds for x
%   xupper  : (n,)-dim vector of upper bounds for x
%   clower  : (m,)-dim vector of lower bounds for x
%   cupper  : (m,)-dim vector of lower bounds for x
%   x0      : (n,)-dim vector of inital point for x
% 
%%

