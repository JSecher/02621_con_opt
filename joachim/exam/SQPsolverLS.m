function [x, lambda, output] = SQPsolverLS(objfun,confun,xlower,xupper,clower,cupper,x0)
% SQPsolverLS SQP solvers using dampend BFGS approximation and line seach.
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
tol_c1 = 0.1;   % Tolerance for line search
non_monotone = true;

% Allocate storage
n = length(x0);
x = reshape(x0, n, 1);      % Make sure x is a collumn vector
[f,df] = objfun(x);         % Compute obj and con function for x0
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
d = [xlower; -xupper; clower; -cupper];
mu = 0;

% Create outputs struct
output.iterations = 0;
output.converged = false;
output.xk = x;
output.stepLengths = zeros(0,1);
output.time_qp = 0;
output.function_calls = 2;

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
    
    %% Do line search for step length  
    % Set starting alpha at max
    alpha = 1.0;
    % Compute reference merit
    c_alpha = -1*confun(x); 
    output.function_calls = output.function_calls + 2;
    c_alpha = [x; -x; c_alpha; -c_alpha]-d;
    absz = abs(z);
    mu = max(absz, 0.5*(mu+absz));
    phi0 = phi(f,mu,c_alpha);
    dphi0 = dphi(df,Deltax,mu,c_alpha); 
    % Run til convergence
    searchingForAlpha = true;
    while searchingForAlpha
        % Compute step canditate
        x_alpha = x + alpha*Deltax;
        % Compute values for step candidate
        f_alpha = objfun(x_alpha);       
        c_alpha = -1*confun(x_alpha);
        output.function_calls = output.function_calls + 2;
        c_alpha = [x_alpha; -x_alpha; c_alpha; -c_alpha]-d;
        % Compute merit of step
        phi_alpha = phi(f_alpha,mu,c_alpha);

        if phi_alpha <= phi0 + tol_c1*dphi0*alpha
            % Accept alpha value
            searchingForAlpha = false;
        else
            % Update alpha value
            alpha_tmp = (phi_alpha-(phi0+alpha*dphi0))/(alpha*alpha);
            alpha_min = -dphi0/(2*alpha_tmp);
            alpha = min(0.9*alpha, max(alpha_min, 0.1*alpha));
        end
    end

    % Use non monotone update strategy if alpha is very small, and it is
    % set
    if all(alpha*Deltax<epsilon) && non_monotone
        alpha = 1.0;
    end

    % Save used step length
    output.stepLengths(output.iterations) = alpha;

    % Update the current point using alpha
    z = z + alpha*(zhat-z);
    x = x + alpha*Deltax;
    
    % Save step
    output.xk(:, output.iterations+1) = x;
    
    % For the quasi Newton update  
    dL = df - z(lid)-z(uid)-dc*z(clid)-dc*z(cuid);

    % Compute function values for next iteration
    [f,df] = objfun(x);
    [c,~,dc,~] = confun(x);
    c = -1*c;
    dc = -1*dc;
    output.function_calls = output.function_calls + 2;
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

% end of function
end

% Below is function definitions

% Define merit function
function [res] = phi(f, mu, c)
    res = f + mu'*abs(min(0,c));
end
% Define derivative of merit function
function [res] = dphi(df,dx,mu,c)
    res = df'*dx-mu'*abs(min(0,c));
end
