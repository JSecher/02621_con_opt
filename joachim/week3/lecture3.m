
%% Problem 5 - Generate random convex quadratic programs (consider how this can be done) and test you program.

N = 10; % Number of tests
maxcoef = 10;
nvar = 4;
ncon = 3;
all_true = true;

for i=1:N
    [H, g, A, b, xtarget] = GeneratecConvexQuadraticProblems(maxcoef, nvar, ncon);
    [xest,lambdaest] = EqualityQPSolver(H,g,A,b);
    if norm(xest - xtarget) > 10e-4
        all_true = false;
        disp('The solution is not correct')
    end
end

if all_true
    disp("All solutions true")
end


%% Problem 6 - Write the sensitivity equations for the equality constrained convex QP

