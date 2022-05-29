function [x, lambda, output] = SQPsolverTR(objfun,confun,xlower,xupper,clower,cupper,x0)
% SQPsolverBFGS SQP solvers using dampend BFGS approximation and Trust Region.
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
epsilon = 1e-4;

% Allocate storage
n = length(x0);
x = reshape(x0, n, 1);      % Make sure x is a collumn vector
[f,df] = objfun(x);         % Compute obj and con function for x0
[c,~,dc,~] = confun(x);
c = -1*c;     % Fix sign
dc = -1*dc;   % Fix sign
m = size(c,1);
B = eye(n);
z = ones(2*m+2*n,1);
lid = 1:m;
uid = (m+1):(2*m);
clid = (2*m+1):(2*m+n);
cuid = (2*m+n+1):(2*(n+m));

% Trust region specific stuff
tr0 = 0.5;      % Initial trust region size
tr = tr0;       % Iterative trust region size
mu_val = 1000;    % Initial penalty value
mu = mu_val * ones(2*m+2*n,1);  % Make a vector for stuff
dL2 = Inf;

% Allocate storange for penalty program
penaltyH = zeros(3*n+2*m);
nm2 = 2*n+2*m;

% Create outputs struct
output.iterations = 0;
output.converged = false;
output.xk = x;
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

    % Define and solve program for with penalty and tr constraints
    penaltyH(1:n,1:n) = 0;
    penaltyH(1:n,1:n) = B;
    penaltyC = [eye(n), -eye(n), dc, -dc, zeros(n,nm2), eye(n), -eye(n); ...
                     eye(nm2), eye(nm2), zeros(nm2,n*2)]';
    penaltyg = [df; mu];
    penaltyd = [xlowerk; -xupperk; ...
                   clowerk; -cupperk; ...
                   zeros(nm2,1);  ...
                   -tr*ones(2*m,1)];
    start = cputime;
    [Deltax,~,flag,out,lambda] = quadprog(penaltyH,penaltyg, ...
                                          -penaltyC,-penaltyd, ...
                                          [],[],[],[],[],options);
    output.time_qp = output.time_qp + (cputime-start);
    
    if flag < 0
        error("Quadprog error: %s", out.message)
    end
    
    % Get only the relevant parts
    zhat = lambda.ineqlin(1:nm2);
    Deltax = Deltax(1:n);
    
    % Compute new penalty
    z_inf = norm(z, "inf");
    mu_val = max(0.5*(mu_val+z_inf), z_inf);
    mu = mu_val * ones(nm2, 1);

    % Compute actual / predicted ratio
    [c_full,~,dc_full,~] = confun(x);  % Current con
    c_full = [x; -x; c_full; -c_full] - penaltyd(1:nm2);
    dc_full = [eye(n), -eye(n), dc_full, -dc_full];
    c_full = -1*c_full;     % Fix sign
    dc_full = -1*dc_full;   % Fix sign
    [c_pred_full] = confun(x+Deltax);  % predicted con
    c_pred_full = [x; -x; c_pred_full; -c_pred_full] - penaltyd(1:nm2);
    c_pred_full = -1*c_pred_full;     % Fix sign
    
    % Compute predicted 
    qmu_pred = f + df'*Deltax + 0.5*Deltax'*B*Deltax + mu'*max(0,-(c_full+dc_full'*Deltax));
    qmu = f+ mu'*max(0,-(c_full));
    phi1 = qmu;

    [f_pred,~] = objfun(x+Deltax); % predicted obj
    phi1_pred = f_pred+mu'*max(0,-c_pred_full);
    
    rho = (phi1-phi1_pred)/(qmu-qmu_pred);
    gamma = min(max( (2*rho-1)^3 + 1, 0.25), 2);
    
    output.function_calls = output.function_calls + 3;
    
    % Adjust trust region accordingly
    if rho > 0    
        % If accepted
        % Update the current point
        z = zhat;
        x = x + Deltax;
        
        % Save step
        output.xk = [output.xk, x];
        
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
    
        r = theta*q+(1-theta)*Bp;
        B = B + r*r'/(Deltax'*r) - Bp*Bp'/pBp;

        % Update trust region
        tr = gamma*tr;
    else
        % Not accepted
        % Update trust region
        tr = gamma*norm(Deltax,"inf");
    end
    
    % Check for convergence
    if norm(dL2, 'inf') < epsilon
        output.converged = true;
    end
end

% Warn if not converged
if ~output.converged
    warning("Max number of iterations reached before convergence")

end


