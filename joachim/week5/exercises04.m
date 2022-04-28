%% Lecture 4
%% 1. Use fmincon in Matlab’s Optimization Toolbox
% a. Solve Himmelblaus optimization constrained optimization problem

himobjfun = @(x,p) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;


x0 = [0;0]; % initial point
xl = [-5;-5]; % lower bounds
xu = [5;5]; % upper bounds
% First we pretend that there is no linear equality/inequality constraints
A = zeros(0,2);
b = zeros(0,1);
Aeq = zeros(0,2);
beq = zeros(0,1);
% Parameters (in this case we do not use parameters)
p = [];
% Call fmincon
options = optimoptions( 'fmincon','Display','none','Algorithm','interior-point');

[x,fval,exitflag,output]=fmincon(himobjfun, x0, A, b, Aeq, beq, xl, xu, @confunhimmelblau, options, p);

x,fval,output

xtrue = x;

% c try different algorithms

algos = {'interior-point', 'sqp', 'sqp-legacy', 'active-set'};

for i=1:length(algos)
    curralgo = string(algos(i));
    options = optimoptions( 'fmincon','Display','none','Algorithm',curralgo);

    [x,fval,exitflag,output]=fmincon(himobjfun, x0, A, b, Aeq, beq, xl, xu, @confunhimmelblau, options, p);
    
    disp(sprintf("\nAlgorithm: %s", curralgo));
    disp(sprintf("MNumber of iterations: %d", output.iterations));
    

   
    if x(1) ~= xtrue(1)
        str = sprintf("x1 is different than interior-point, err: %0.5e", x(1)-xtrue(1));
        disp(str)
    end
    if x(2) ~= xtrue(2)
        str = sprintf("x2 is different than interior-point, err: %0.5e", x(2)-xtrue(2));
        disp(str)
    end



end





function [c,ceq] = confunhimmelblau(x, p)
    c = zeros(2,1);     % Compute nonlinear inequalities at x.
    ceq = zeros(0,1);   % Compute nonlinear equalities at x.
    tmp = x(1)+2;
    c(1,1) = -(tmp*tmp - x(2));
    c(2,1) = -(-4*x(1) + 10*x(2));
end