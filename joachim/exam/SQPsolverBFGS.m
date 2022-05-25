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
n = length(x0);
x = reshape(x0, n, 1);      % Make sure x is a collumn vector
[~,df] = objfun(x);         % Compute obj and con function for x0
[c,~,dc,~] = confun(x);
c = -1*c;
dc = -1*dc;
m = size(c,1);
B = eye(n);
z = ones(2*m+2*n,1);
lid = 1:m;
uid = (m+1):(2*m);
clid = (2*m+1):(2*m+n);
cuid = (2*m+n+1):(2*(n+m));


% Create outputs struct
output.iterations = 0;
output.converged = false;
output.xk = x;
output.time_qp = 0;

% Define options for quadprog
options = optimset('Display', 'off');

% Run until convergence or max iter
while (output.iterations < maxiter) && ~output.converged
    
    output.iterations = output.iterations + 1;

    % Set lower and upper bounds iteration k 
    xlowerk = -x + xlower;
    xupperk = -x + xupper;
    clowerk = -c + clower;
    cupperk = -c + cupper;
    
    start = cputime;
    [Deltax,~,flag,out,lambda] = quadprog(B,df, ...
                                     -[dc'; -dc'],-[clowerk;-cupperk], ...
                                     [],[], ...
                                     xlowerk,xupperk,[], ...
                                     options);
    output.time_qp = output.time_qp + (cputime-start);
    
    if flag < 0
        error("Quadprog error: %s", out.message)
    end
    zhat = [lambda.lower; lambda.upper; lambda.ineqlin];
    
    % Update the current point
    z = z + (zhat-z);
    x = x + Deltax;
    
    % Save step
    output.xk(:, output.iterations+1) = x;
    
    % For the quasi Newton update  
    dL = df - z(lid)-z(uid)-dc*z(clid)-dc*z(cuid);

    % Compute function values for next iteration
    [~,df] = objfun(x);
    [c,~,dc,~] = confun(x);
    c = -1*c;
    dc = -1*dc;
    %% Dampend BFGS Update
    
    % Compute Quasi Newton update of the hessian
    dL2 = df - z(lid)-z(uid)-dc*z(clid)-dc*z(cuid);

    
    % Compute values used multiple times for for BFGS update 
    q = dL2-dL;
    Bp = (B*Deltax);
    pBp = Deltax'*Bp;
    
    % Compute appropiate theta value
    if Deltax'*q < 0.2*pBp
        theta = (0.8*pBp)/(pBp-Deltax'*q);
    else
        theta = 1;
    end

    r = theta*q+(1-theta)*(Bp);
    B = B + r*r'/(Deltax'*r) - Bp*Bp'/pBp;
    
    % Check for convergence
    if norm(dL2, 'inf') < epsilon
        output.converged = true;
    end
end

% Warn if not converged
if ~output.converged
    warning("Max number of iterations reached before convergence")

end


