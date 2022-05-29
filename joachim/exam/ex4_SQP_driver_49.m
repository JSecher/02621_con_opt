%% Problem 4.9 Driver
%% Clean up
clear
close all

Himmel_plots = true;
Rosen_plots = true;
rosen_con_case = 1;

% Import CasADi
if isfolder('../../casadi-v3.5.5')
    addpath('../../casadi-v3.5.5')
elseif isfolder('../../../casadi-v3.5.5')
    addpath('../../../casadi-v3.5.5')
else
    error("Can not find casadi")
end
import casadi.*
disp("Casadi import succesfull")


%% Problem 4.9 - Himmelblau
if Himmel_plots
%%
disp("Solving Himmelblaus Test problem all own solvers");

x0s = [[.5; 4.0], [-4.0; 0.5]];


own_solvers = ["bfgs", "line", "trust"];

for j=1:length(x0s)

    fig = figure("Name", "Himmel - solutions for x0 ", 'Position', [100, 100, 800,800]);
    hold on
    %Create contour plot and constraints
    [cfig, conFigs] = contourHimmel(true);
    x0 = x0s(:,j);    % Initial point
    xl = [-5; -5];      % Lower bound for x
    xu = [5; 5];        % Upper bound for x
    cl = [0; 0];        % Lower bound for constraints 
    cu = [54; 70];      % Upper bound for constraints

    legs = [];
    tits = [];
    for i=1:length(own_solvers)
        solv = own_solvers(i);
        fprintf("Using solver: %s\n", solv)
       
        
        tstart = cputime;
        [sol_own, obj_own, lambda, output_own] = SQPsolver(@objfungradHimmelblau, ...
                                                    @confungradHimmelblau1, ...
                                                    xl, xu, ...
                                                    cl, cu, ...
                                                    x0, solv);
                                                               
        time_own = cputime - tstart;
            
        fprintf('\t %s solution: [%.5e, %.5e], objective: %.5e, Iter: %.5e\n', solv, sol_own(1), sol_own(2), obj_own, output_own.iterations);

        colors = ['r', 'g', 'y'];
        names = ["BFGS", "Line Search", "Trust Region"];
        col = colors(i);
        % Add points of interest
        traceIterations(output_own.xk, col, '-');
        ow = plotPoint(sol_own(1),sol_own(2), "sol", col, 18);
        legs = [legs, ow];
        tits = [tits, names(i)];
    end
        
    % Solve problem using fmincon
    options = optimoptions('fmincon',... 
                           'Display','none',... 
                           'Algorithm','interior-point');

    % For fmincon
    tstart = cputime;
    [sol_fmin,fval,exitflag,output] = fmincon( ...
                                            @objfunHimmelblau, x0, ...
                                            [], [], [], [], ...
                                            xl, xu, ...
                                            @confunHimmelblau, ...
                                            options);
    time_fmincon = cputime - tstart;

    % Print results
    fprintf('Found solutions for x0 = [%.2f, %.2f] ::\n', x0(1), x0(2)); 
    fprintf('\t fmincon solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_fmin(1), sol_fmin(2), fval, time_fmincon);
    
    
    %data(j,:) = [fval, sol_fmin', time_fmincon, nan, ...
    %             obj_cas, sol_cas', time_cas, nan, mean(sqrt((sol_fmin-sol_cas).^2))];

        
    % Add points of interest
    %fm = traceIterations([x0, sol_fmin], "r", '-');
    
    intp = plotPoint(x0(1),x0(2), "int");
    solp_fm = plotPoint(sol_fmin(1),sol_fmin(2), "sol", "b", 10);
    legs = [legs, intp, solp_fm];
    tits = [tits, "Initial point", "Fmincon"];

    legend(legs,tits, "Location","southeast")
    hold off
    %savefigpdf(fig, sprintf('ex4_9_solve_x0_%d_himmel', j), 4);
end




end


%% Problem 4.9 - Solve Rosenbrock 
if Rosen_plots
%%
disp("Solving Rosenbrock Test problem all own solvers");

% Select constraint for rosenbrok, 1 = circle, 2 = box
switch rosen_con_case
    case 1
        con_type_name = "Circle";
        x0s = [[0.0; 0.0], [-0.5; 0], [-0.5; -0.5]];
    case 2
        con_type_name = "Box";
        x0s = [[0.5; 1.0], [-1.25; 0.5], [-0.5; 1], [0.75; 0.6]];
    otherwise
        error("Unknown case for rosenbrock")
end

own_solvers = ["bfgs", "line", "trust"];

for j=1:length(x0s)

    fig = figure("Name", sprintf("Rosenbrock - solutions for x0 %s constraints", con_type_name), 'Position', [100, 100, 800,800]);
    hold on
    %Create contour plot and constraints
    [cfig, conFigs] = contourRosen(rosen_con_case);

    legs = [];
    tits = [];
    for i=1:length(own_solvers)
        solv = own_solvers(i);
        fprintf("Using solver: %s\n", solv)
        switch rosen_con_case
            case 1
                x0 = x0s(:,j);    % Initial point
                xl = [-1; -1];      % Lower bound for x
                xu = [1; 1];        % Upper bound for x
                cu = [1;100];
                cl = [0;-100];
                % For fmincon
                tstart = cputime;
                [sol_own,obj_own,~,output_own] = SQPsolver(@objfungradRosenbrock, ...
                                                            @confungradRosenbrock,...
                                                            xl, xu, ...
                                                            cl, cu,...
                                                            x0, solv);
                                                           
                time_own= cputime - tstart;
            case 2
                x0 = x0s(:,j);    % Initial point
                xl = [-1.2; 0.5];      % Lower bound for x
                xu = [1.2; 1.5];        % Upper bound for x
                cl = [-100; -100];
                cu = [100; 100];
                % fmincon
                tstart = cputime;
                [sol_own,obj_own,~,output_own] = SQPsolver(@objfungradRosenbrock, ...
                                                            @confungradRosenbrock_part2,...
                                                            xl, xu, ...
                                                            cl, cu,...
                                                            x0, solv);
                                                           
                time_own = cputime - tstart;
            otherwise
                error("Unknown case for rosenbrock")
        end
        fprintf('\t %s solution: [%.5e, %.5e], objective: %.5e, Iter: %.5e\n', solv, sol_own(1), sol_own(2), obj_own, output_own.iterations);

        colors = ['r', 'g', 'y'];
        names = ["BFGS", "Line Search", "Trust Region"];
        col = colors(i);
        % Add points of interest
        traceIterations(output_own.xk, col, '-');
        ow = plotPoint(sol_own(1),sol_own(2), "sol", col, 18);
        legs = [legs, ow];
        tits = [tits, names(i)];
    end
        
    % Solve problem using fmincon
    options = optimoptions('fmincon',... 
                           'Display','none',... 
                           'Algorithm','interior-point');

    switch rosen_con_case
        case 1
            x0 = x0s(:,j);    % Initial point
            xl = [-1; -1];      % Lower bound for x
            xu = [1; 1];        % Upper bound for x
            % For fmincon
            tstart = cputime;
            [sol_fmin,fval,exitflag,output] = fmincon(@objfunRosenbrock, ...
                                                       x0, ...
                                                       [], [], [], [], ...
                                                       xl, xu, ...
                                                       @confunRosenbrock, ...
                                                       options);
            time_fmincon = cputime - tstart;
            
        case 2
            x0 = x0s(:,j);    % Initial point
            xl = [-1.5; 0.5];      % Lower bound for x
            xu = [1.5; 1.5];        % Upper bound for x
            % fmincon
            tstart = cputime;
            [sol_fmin,fval,exitflag,output] = fmincon(@objfunRosenbrock, ...
                                                       x0, ...
                                                       [], [], [], [], ...
                                                       xl, xu, ...
                                                       @confunRosenbrock_part2, ...
                                                       options);
            time_fmincon = cputime - tstart;
            
        otherwise
            error("Unknown case for rosenbrock")
    end


    % Print results
    fprintf('Found solutions for x0 = [%.2f, %.2f] ::\n', x0(1), x0(2)); 
    fprintf('\t fmincon solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_fmin(1), sol_fmin(2), fval, time_fmincon);
    
    
    %data(j,:) = [fval, sol_fmin', time_fmincon, nan, ...
    %             obj_cas, sol_cas', time_cas, nan, mean(sqrt((sol_fmin-sol_cas).^2))];

        
    % Add points of interest
    %fm = traceIterations([x0, sol_fmin], "r", '-');
    
    intp = plotPoint(x0(1),x0(2), "int");
    solp_fm = plotPoint(sol_fmin(1),sol_fmin(2), "sol", "b", 10);
    legs = [legs, intp, solp_fm];
    tits = [tits, "Initial point", "Fmincon"];

    legend(legs,tits, "Location","southwest")
    hold off
    savefigpdf(fig, sprintf('ex4_9_solve_x0_%d_test_rosen_%s', j, lower(con_type_name)), 4);
end




end

%% Function definition 

function f = objfunHimmelblau(x,p)
    % Function giving the objective of the Himmelblau problem
    % Function taken from: Slide 8, Lecture 01B 
    tmp1 = x(1)*x(1)+x(2)-11;
    tmp2 = x(1)+x(2)*x(2)-7;
    f = tmp1*tmp1 + tmp2*tmp2;
end

function [c,ceq] = confunHimmelblau(x,p)
    % Function giving the constraints for the Himmelblau problem
    % Function taken from: Slide 8, Lecture 01B 
    c = zeros(2,1);
    ceq = zeros(0,1);

    % Inequality constraints c(x) <= 0
    tmp = x(1)+2;
    c(1,1) = -(tmp*tmp - x(2));
    c(2,1) = -(-4*x(1) + 10*x(2));
end

function [f,dfdx] = objfungradHimmelblau(x,p)
    % Function giving the objective and gradients of the Himmelblau problem
    % Function taken from: Slide 8, Lecture 01B 
    tmp1 = x(1)*x(1)+x(2)-11; 
    tmp2 = x(1)+x(2)*x(2)-7;
    f = tmp1*tmp1 + tmp2*tmp2;
    % compute the gradient of f
    if nargout > 1
        dfdx = zeros(2,1);
        dfdx(1,1) = 4*tmp1*x(1) + 2*tmp2; 
        dfdx(2,1) = 2*tmp1 + 4*tmp2*x(2);
    end
end

function [c,ceq,dcdx,dceqdx] = confungradHimmelblau1(x,p) 
    % Function giving the constraints and gradients of the Himmelblau problem
    % Function taken from: Slide 8, Lecture 01B 
    c = zeros(2,1);
    ceq = zeros(0,1);
    % Inequality constraints c(x) <= 0 
    tmp = x(1)+2;
    c(1,1) = -(tmp*tmp - x(2));
    c(2,1) = -(-4*x(1) + 10*x(2));
    % Compute constraint gradients
    if nargout > 2
        dcdx = zeros(2,2); 
        dceqdx = zeros(2,0);
        dcdx(1,1) = -2*tmp; % dc1dx1
        dcdx(2,1) = 1.0; % dc1dx2
        dcdx(1,2) = 4; % dc1dx1
        dcdx(2,2) = -10; % dc1dx2
    end
end

function f = objfunRosenbrock(x,p)
    % Function giving the objective of the Rosenbroc problem,
    % objective function
    tmp1 = (x(2)-x(1)*x(1));
    tmp2 = tmp1*tmp1;
    f = 100 * tmp2 + (1-x(1))*(1-x(1));
end

function [f,dfdx] = objfungradRosenbrock(x,p)
    % Function giving the constraints and gradients of Rosenbroc problem,
    % objective function
    tmp1 = (x(2)-x(1)*x(1));
    tmp2 = tmp1*tmp1;
    f = 100 * tmp2 + (1-x(1))*(1-x(1));
    % compute the gradient of f
    if nargout > 1
        dfdx = zeros(2,1);
        dfdx(1,1) = -400*(-x(1)^2 + x(2))*x(1) - 2 + 2*x(1);
        dfdx(2,1) = -200*x(1)^2 + 200*x(2);
    end
end

function [c,ceq] = confunRosenbrock(x,p)
    % Function giving the constraints of Rosenbrock problem,
    % constraints 1, Unit circle
    c = zeros(1,1);
    ceq = zeros(0,1);

    % Inequality constraints c(x) <= 0
    c(1,1) = x(1)^2 + x(2)^2 - 1;
end

function [c,ceq] = confunRosenbrock_part2(x,p)
    % Function giving the constraints of Rosenbrock problem,
    % constraints 2, Box
    c = zeros(0,1);
    ceq = zeros(0,1);
end

function [c,ceq,dcdx,dceqdx] = confungradRosenbrock(x,p) 
    % Function giving the constraints and gradients of Rosenbroc problem,
    % constraints 1, unit circle
    c = zeros(2,1);
    ceq = zeros(0,1);
    % Inequality constraints c(x) <= 0 
    c(1,1) = x(1)^2 + x(2)^2 - 1;
    c(2,1) = x(1) + x(2);
    % Compute constraint gradients
    if nargout > 2
        dcdx = zeros(2,2); 
        dceqdx = zeros(2,0);
        dcdx(1,1) = 2*x(1); % dc1dx1
        dcdx(2,1) = 2*x(2); % dc1dx2
        dcdx(1,2) = 1; % dc1dx1
        dcdx(2,2) = 1; % dc1dx2
    end
end

function [c,ceq,dcdx,dceqdx] = confungradRosenbrock_part2(x,p) 
    % Function giving the constraints and gradients of Rosenbroc problem,
    % constraints 2, Box
    c = zeros(0,1);
    ceq = zeros(0,1);

    % added fake
    c = zeros(2,1);
    ceq = zeros(0,1);
    c(1,1) = x(1) + x(2);
    c(2,1) = x(1) + x(2);
    % Inequality constraints c(x) <= 0 
    % Compute constraint gradients
    if nargout > 2
        dcdx = zeros(2,2); 
        dceqdx = zeros(2,0);
        dcdx(1,1) = 1; % dc1dx1
        dcdx(2,1) = 1; % dc1dx2
        dcdx(1,2) = 1; % dc1dx1
        dcdx(2,2) = 1; % dc1dx2
    end
end


function [c,ceq,dcdx,dceqdx] = confungradRosenbrock_part3(x,p) 
    % Function giving the constraints and gradients of Rosenbroc problem,
    % constraints 2, Box
    c = zeros(0,1);
    ceq = zeros(0,1);

    % added fake
    c = zeros(2,1);
    ceq = zeros(0,1);
    c(2,1) = 2*(x(1)*x(1)) - x(2);
    c(1,1) = x(1) + x(2);
    % Inequality constraints c(x) <= 0 
    % Compute constraint gradients
    if nargout > 2
        dcdx = zeros(2,2); 
        dceqdx = zeros(2,0);
        dcdx(1,1) = 4*x(1); % dc1dx1
        dcdx(2,1) = 4*x(2); % dc1dx2
        dcdx(1,2) = 1; % dc1dx1
        dcdx(2,2) = 1; % dc1dx2
    end
end

function [h] = traceIterations(xks, color, linetype)
% Display the iterations for a SQP algortihm for two variables
    if nargin<2
        color = 'r';
    end
    if nargin<3
        linetype = '-.';
    end
    spec = sprintf("%so%s", linetype, color);
    Nvals = size(xks, 2);
    plot(xks(1,[1, Nvals]),xks(2,[1,Nvals]),"MarkerFaceColor", color,'markersize',8,'LineStyle', 'none' );
    h = plot(xks(1,:),xks(2,:),spec,'linewidth',2,'MarkerSize',5);
end