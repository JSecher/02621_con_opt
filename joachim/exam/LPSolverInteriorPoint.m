function [x,y,z,s,iter] = LPSolverInteriorPoint(g,A,b,l,u,x0,y0,z0,s0,tol,maxiter, eta)
% LPSolverInteriorPoint Mehrotra s predictor-corrector primal-dual
% interior-point algorithm LP with box constraints
%
%   Detailed explanation goes here
%
% Inputs:
%   g       : a (n,) dimensional vector defining the linear part of the objective function
%   A       : a (n,m) dimensional matrix giving the rhs of the equality constrains
%   b       : a (m,) dimensional vector defining the rhs of the equality constrains
%   l       : a (n,) dimensional vector defining lower bounds of x
%   u       : a (n,) dimensional vector defining upper bounds of x
%   x0      : a (n,) dimensional vector with initial guess for variables
%   y0      : a (m,) dimensional vector with initial guess for equality lagrange multipliers
%   z0      : a (2n,) dimensional vector with initial guess for inequality
%               lagrange multipliers. All values must be > 0
%   s0      : a (n,) dimensional vector with initial guess for the slack
%               variables. All values must be > 0
%   tol     : a float for tolerance for stopping criteria
%   maxiter : an integer giving the max number of iterations allowed
%
% Outputs:
%   x           : the solution as a (n,) dimensional vector with 
%   y           : a (m,) dimensional vector with equality lagrange multipliers
%   z           : a (2n,) dimensional vector with inequality lagrange multipliers
%   s           : a (n,) dimensional vector with the slack variables
%   iterations  : a integer giving the number of iterations until convergens
%% 

% Constants 
n = length(x0);
m = length(y0);
m_c = length(z0);

if nargin < 10
    tol = 1e-9;
end
if nargin < 11
    maxiter = 1000;
end
if nargin < 12
    eta = 0.995;
end

% Avoid singular matrix operations 
while (any(s0==0))
    x0 = x0 + 1e-6;
    s0 = [x0-l; -x0+u];
end

debug = true;
if debug; disp("Debug ON"); end

%% --- Improve the starting point using Heuristic ----
% Compute Initial residuals 
[rL,rA,rC,rSZ] = ComputeResiduals(g,A,b,l,u,x0,y0,z0,s0);

% Compute the affine step from the newton direction
%rCs = rC-s0;
rCmod = s0;
[~, ~, dzAff, dsAff, ~] = NewtonDirection(rA,rL,rC,rCmod,A,z0,s0,[]);

% Simply reuse x and y
x = x0;
y = y0;
% Set z and s from the Affine step
z = max(1, abs(z0+dzAff));
s = max(1, abs(s0+dsAff));


%% ---- Get ready for first iteration from updated starting guess ----
% Update Initial residuals 
[rL,rA,rC,rSZ] = ComputeResiduals(g,A,b,l,u,x,y,z,s);
    
% Initial dual gap used for convergence check
mu0 = (z'*s)/m_c;
mu = mu0;

% Run Interior-point until convergence
converged = false;
iter = 0;
while (iter < maxiter && ~converged)
    iter = iter + 1;
    if debug
        fprintf("Iter %d, mu = %e\n", iter, mu)
    end
    %% ---- Predictor step ----
    % Compute the affine direction from the newton direction 
    %rCs = rC-s;
    rCmod = s;
    [dxAff, dyAff, dzAff, dsAff, R] = NewtonDirection(rA,rL,rC,rCmod,A,z,s,[]);

    % Get affine-scaling step length
    alphaAff = ComputeStepSize(z, s, dzAff, dsAff);
    
    % Compute the a ffine duality gap
    muAff = (z+alphaAff*dzAff)'*(s+alphaAff*dsAff)/m_c;
    
    % Compute the centering parameter
    sigma = (muAff/mu)^3;
     
    %% ---- Corrector step ----
    
    % Set up right hand sides
    %rCs = rC - (s + dzAff.*dsAff./z - sigma*mu*ones(m_c,1)./z);
    rCmod = (s + dzAff.*dsAff./z - sigma*mu*ones(m_c,1)./z);
    % Get the corrector direction
    [dxCC, dyCC, dzCC, dsCC, ~] = NewtonDirection(rA,rL,rC,rCmod,A,z,s,R);
    
    %% ---- Compute step size and take step ----
    dx = dxCC;
    dy = dyCC;
    dz = dzCC;
    ds = dsCC;
    
    alpha = ComputeStepSize(z, s, dz, ds);
    
    % Update of position 
    alphaBar = eta * alpha;
    x = x + alphaBar * dx;
    y = y + alphaBar * dy;
    z = z + alphaBar * dz;
    s = s + alphaBar * ds;

    %% ---- Update for next iter and check for convergence ----
    % Update residuals and compute the dual gap
    [rL,rA,rC,rSZ] = ComputeResiduals(g,A,b,l,u,x,y,z,s);
    mu = (z'*s)/m_c;

    % Check convergence
    converged = ConvergenceCheck(tol, mu, mu0);
end

if iter >= maxiter
    warning("QPSolverMehrotraInteriorPoint did not converge within %d iterations", maxiter);
end

end


function [dx, dy, dz, ds, R] = NewtonDirection(rA,rL,rC,rCmod,A,z,s,R)
% NewtonDirection Compute the affine Newton Direction in the Interior-point algorihtm
%
% Inputs:
%   x0      : a with initial guess for variables
%   y0      : a with initial guess for equality lagrange multipliers
%   z0      : a with initial guess for inequality
%   s0      : a vecotr with initial guess for the slack variables
%
% Outputs:
%   dx    : the newton direction for x
%   dz    : the newton direction for z
%   ds    : the newton direction for s
%   R     : R matrix from Cholesky factorization 
%% 

[n,m] = size(A);

zs = z./s;
diagHinv = 1./(zs(1:n)+zs(n+1:end));

rCs = rC - rCmod;
rLbar = rL - zs(1:n).*rCs(1:n) + zs(n+1:end).*rCs(n+1:end);

% Compute Cholesky-factorization if R is not provided
if isempty(R)
    nf = A' * (diagHinv .*  A);
    R = chol(nf);
end

dy = R \ (R' \ (rA + A' * (diagHinv .* rLbar)));
dx = diagHinv .* (-rLbar + A*dy);
dz = -zs.*[dx; -dx] + zs.*rCs;
ds = -rCmod-(s./z).*dz;

end

function [rL,rA,rC,rSZ] = ComputeResiduals(g,A,b,l,u,x,y,z,s)
% StartingPointHuristic Heuristic for improving the starting point of
% of a interior-point algorithm 
% Computed from the formular:
%      rL = g - Ay - Cz 
%      rA = -A'x + b
%      rC = -C'x + s + d \\
%      rSZ = SZe
%
% Outputs:
%   rL    : Lagrangian gradient
%   rA    : Equality Constraint
%   rC    : Inequality Constraint
%   rSZ   : Complementarity
%% 

rL = g-A*y-(z(1:length(x))-z(length(x)+1:end));
rA = b-A'*x;
rC = s + [l; -u] - [x; -x];
rSZ = s.*z;

end

function [alpha] = ComputeStepSize(z, s, dz, ds)
% ComputeStepSize Compute the step size in the interior point alogrthm 
% Compute the largest alpha such that: z + alpha*∆zaff  >=  0 and s + alpha*∆saff >=  0

% Find largest alpha for z
idx = find(dz < 0.0);
alphaz = min([1.0; -z(idx,1)./dz(idx,1)]);
% Find largest alpha for s
idx = find(ds < 0.0);
alphas = min([1.0; -s(idx,1)./ds(idx,1)]);
% Compute the largest alpha such that: z + alpha*∆zaff  >=  0 and s + alpha*∆saff >=  0
alpha = min(alphaz,alphas);

end

function [converged] = ConvergenceCheck(tol,mu,mu0)
% ConvergenceCheck Checks stopping criteria of interior-point algorithm 
% 
% See Gertz and Wright (2001)
%
% Inputs:
% 
% Outputs:
%   converged   : a bool, true if convergence has been reached
%% 

% Converged, only mu is implemented
converged = (mu <= tol*1e-2*mu0);

end