function [x, lambda, output] = SQPsolverBFGS(objfun,confun,xlower,xupper,clower,cupper,x0)
% SQPsolverBFGS SQP solvers using dampend BFGS approximation.
%
% Solves NLP's on the form  
%
%   \min_{x} f(x) 
%    s.t.    cl <= c(x) <= cu  
%            xl <= x <= xu
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
%
% Outputs:
%   x       : a (n,) dimensional vector of found solution
%   lambda  : a (m,) dimensional vector of lagrange multipliers
%   output  : a object with performance and iteration information
%
%%

% Define constants 
maxiter = 100;
epsilon = 1e-9;

% Allocate storage
x = x0;
[~,df] = objfun(x);
[c,dc] = confun(x);
n = length(x);
m = size(c,1);
B = eye(n);

% Define options for quadprog
options = optimset('Display', 'off');

% Run until convergence or max iter
iter = 0;
converged = false;
while (iter < maxiter) && ~converged
    
    % Update lower and upper bounds for the quadrastart = cputime; approximation
    lk = -x+l;
    uk = -x+u;
    clk = -c+cl;
    cuk = -c+cu;
    
    start = cputime;
    [Deltax,~,~,~,lambda] = quadprog(B,df,-[dc'; -dc'],-[clk;-cuk],[],[],lk,uk,[], options);
    time = cputime-start;

    zhat = [lambda.lower; lambda.upper; lambda.ineqlin];
    
    pz = zhat-z;

    % Update the current point
    z = z + pz;
    x = x + Deltax;
    
    % For the quasi Newton update  
    dL = df - (z(lid)-z(uid)+dc*z(clid)-dc*z(cuid));

    % Compute function values for next iteration
    [~,df] = feval(obj,x);
    [c,dc] = feval(con,x);
    
    %% Dampend BFGS Update
    
    % Compute Quasi Newton update of the hessian
    dL2 = df - (z(lid)-z(uid)+dc*z(clid)-dc*z(cuid));
    
    % Compute set values
    Deltax = pk;
    q = dL2-dL;
    Bp = (B*Deltax);
    pBp = Deltax'*Bp;
    
    if Deltax'*q < 0.2*pBp
        theta = (0.8*pBp)/(pBp-Deltax'*q);
    else
        theta = 1;
    end

    r = theta*q+(1-theta)*(Bp);
    B = B + r*r'/(Deltax'*r) - Bp*Bp'/pBp;
    
    % Check for convergence
    if norm(dL2, 'inf') < epsilon
        converged = true;
    end
end

if ~converged
    warning("Max number of iterations reached before convergence")

end


