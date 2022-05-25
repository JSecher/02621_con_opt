%% Problem 4 Driver
%% Clean up
clear
close all

runContourPlot_44 = false;
runSolveTest_45= false;
runBFGS_46 = true;

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


%% Problem 4.4 - Plot Himmelblau
if runContourPlot_44
%%
disp("Generating and saving contour plot");
% Define points of interest
min_points = [ 3.0, 2.0; ...
              -0.2983476136,2.895620844; ...
              -3.654605171,2.737718273; ...
              -3.548537388,-1.419414955;];

max_points = [-0.4869360343, -0.1947744137; ...
               3.216440661,  1.286576264; ...
              -1.424243078,  0.3314960331];

saddle_points = [0.08667750456, 2.884254701; ...
                 -3.073025751, -0.08135304429];

fig = figure("Name", "Himmelblau - POI", 'Position', [100, 100, 800,800]);
hold on
% Create contour plot and constraints
[cfig, conFigs] = contourHimmel(true);

% Add points of interest
mx = plotPoint(max_points(:,1),max_points(:,2), "max");
mn = plotPoint(min_points(:,1),min_points(:,2), "min");
sd = plotPoint(saddle_points(:,1),saddle_points(:,2), "sad");
if isempty(conFigs)
    legend([mx,mn,sd],{'Local Maximum', 'Local Minimum', 'Saddle Point'})
else
    legend([mx,mn,sd, conFigs(1)],{'Local Maximum', 'Local Minimum', 'Saddle Point', 'Infeasible region'})
end

hold off
savefigpdf(fig, "ex4_himmelblau_poi", 4);

end
%% Problem 4.5 - Solve Himmelblau
if runSolveTest_45
%%
disp("Solving Test problem using fmincon and Casadi");

x0s = [[0.0; 0.0], [-4.0; 0], [-4; 1]];

for j=1:length(x0s)
    %sympref('FloatingPointOutput',1);
    x0 = x0s(:,j);    % Initial point
    xl = [-5; -5];      % Lower bound for x
    xu = [5; 5];        % Upper bound for x
    cl = [0; 0];        % Lower bound for constraints 
    cu = [54; 70];      % Upper bound for constraints
    
    % Solve problem using fmincon
    options = optimoptions('fmincon',... 
                           'Display','none',... 
                           'Algorithm','interior-point');
    tstart = cputime;
    [sol_fmin,fval,exitflag,output] = fmincon(@objfunHimmelblau, ...
                                               x0, ...
                                               [], [], [], [], ...
                                               xl, xu, ...
                                               @confunHimmelblau, ...
                                               options);
    time_fmincon = cputime - tstart;
    
     % Call fmincon
    options = optimoptions('fmincon',... 
                           'SpecifyObjectiveGradient',true,... 
                           'SpecifyConstraintGradient',true,... 
                           'Display','none',... 
                           'Algorithm','interior-point');
    tstart = cputime;
    [sol_fmin_grad,fval_grad,exitflag_grad,output_grad]=fmincon( ...
                                            @objfungradHimmelblau, x0, ...
                                            [], [], [], [], ...
                                            xl, xu, ...
                                            @confungradHimmelblau1, ...
                                            options);
    time_fmincon_grad = cputime - tstart;
                                    
    % Define problem for casadi
    x1 = SX.sym('x1');  % Define variables
    x2 = SX.sym('x2');
    % Define problem
    hbprob = struct('x',[x1; x2], ...
                    'f',(x1^2+x2-11)^2+(x1+x2^2-7)^2, ...
                    'g',[(x1+2)^2-x2; -4*x1+10*x2] ...
                    );
    
    % Solve the problem with casadi
    options = struct;
    options.ipopt.print_level = 0;
    options.print_time = 0;
    tstart = cputime;
    S = nlpsol('S', 'ipopt', hbprob, options);
    r = S('x0', x0, 'lbg',0,'ubg',inf);
    time_cas = cputime - tstart;
    sol_cas = full(r.x);
    obj_cas = full(r.f);
    
    % Print results
    fprintf('Found solutions for x0 = [%.2f, %.2f] ::\n', x0(1), x0(2)); 
    fprintf('\t fmincon solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_fmin(1), sol_fmin(2), fval, time_fmincon);
    fprintf('\t fmincon grad solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_fmin_grad(1), sol_fmin_grad(2), fval_grad, time_fmincon_grad);
    fprintf('\t Casadi solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_cas(1), sol_cas(2), obj_cas, time_cas);
    
    fig = figure("Name", sprintf("Himmelblau - Solution for x0=[%+.1f, %+.1f]",x0(1),x0(2)), 'Position', [150, 150, 600, 600]);
    hold on
    % Create contour plot and constraints
    [cfig, conFigs] = contourHimmel(true);
    
    % Add points of interest
    intp = plotPoint(x0(1),x0(2), "int");
    solp = plotPoint(sol_fmin(1),sol_fmin(2), "sol", 'y', 18);
    solp_grad = plotPoint(sol_fmin_grad(1),sol_fmin_grad(2), "min", 'r', 10);
    solp_cas = plotPoint(sol_cas(1),sol_cas(2), "gen", 'g', 8);
    legend([intp,solp,solp_grad, solp_cas],{'Initial Point', 'fmincon solution', 'fmincon gradient solution','CasADi solution'})
    hold off
    savefigpdf(fig, sprintf("ex4_5_himmelblau_x0=%+.0f_%+.0f",x0(1),x0(2)), 4);

    data(j,:) = [fval, sol_fmin', time_fmincon, nan...
                 fval_grad, sol_fmin_grad', time_fmincon_grad, nan ...
                 obj_cas, sol_cas', time_cas];

end

input.data = data';
% Set column labels (use empty string for no label):
input.tableColLabels = sprintfc("$x_0=[%.1f, %.1f]$", x0s')';
% Set row labels (use empty string for no label):
input.tableRowLabels = {'$f(x)_{\textit{fmincom}} =$', ...
                        '$x_{\textit{fmincom}} = $', '', ...
                        '$\text{time}_{\textit{fmincom}}\, [s] =$', ...
                        '', ...
                        '$f(x)_{\textit{fmincom grad}} =$', ...
                        '$x_{\textit{fmincom  grad}} =$', '', ...
                        '$\text{time}_{\textit{fmincom  grad}}\, [s] =$', ...
                        '', ...
                        '$f(x)_{\textit{CasADi}} =$', ...
                        '$x_{\textit{CasADi}} = $', '', ...
                        '$\text{time}_{\textit{CasADi}}\, [s] =$'};
% Set the row format of the data values 
%input.dataFormatMode = 'row';
input.dataFormat = {'%.5f'};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'r';
% Switch table borders on/off:
input.booktabs = 1;
% LaTex table caption:
input.tableCaption = sprintf('Found solution for different inital points for the Himmelblau test problem.');
% LaTex table label:
input.tableLabel = 'ex4_solve_test_himmel';
input.makeCompleteLatexDocument = 0;
input.dataNanString = '';
input.tablePlacement = '!ht';
% Now call the function to generate LaTex code:
latex = latexTable(input);
savelatexTable(latex, input.tableLabel, 4);


end
%% Problem 4.6 - Dampend BFGS
if runBFGS_46 

clear data
x0s = [[0.0; 0.0],[1.0; 2.0], [-4.0; 0], [-4; 1]];

for j=1:length(x0s)
    %sympref('FloatingPointOutput',1);
    x0 = x0s(:,j);    % Initial point
    
    xl = [-5; -5];      % Lower bound for x
    xu = [5; 5];        % Upper bound for x
    cl = [0; 0];        % Lower bound for constraints 
    cu = [47; 70];      % Upper bound for constraints
    
    tstart = cputime;
    [sol_bfgs, obj, lambda, output] = SQPsolver(@objfungradHimmelblau, ...
                                                @confungradHimmelblau1, ...
                                                xl, xu, ...
                                                cl, cu, ...
                                                x0, 'bfgs');
    t_bfgs_total = cputime - tstart;
    
    
    % Compare with fmin con
    options = optimoptions('fmincon',... 
                           'SpecifyObjectiveGradient',true,... 
                           'SpecifyConstraintGradient',true,... 
                           'Display','none',... 
                           'Algorithm','interior-point');
    tstart = cputime;
    [sol_fmin_grad,fval_grad,exitflag_grad,output_grad]=fmincon( ...
                                            @objfungradHimmelblau, x0, ...
                                            [], [], [], [], ...
                                            xl, xu, ...
                                            @confungradHimmelblau1, ...
                                            options);
    time_fmincon_grad = cputime - tstart;
    
    fprintf('Found solutions for x0 = [%.2f, %.2f] ::\n', x0(1), x0(2)); 
    fprintf('\t BFGS solution: [%.5e, %.5e], objective: %.5e, time: %.5e, iter: %d\n', sol_bfgs(1), sol_bfgs(2), obj, t_bfgs_total, output.iterations);
    fprintf('\t fmincon grad solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_fmin_grad(1), sol_fmin_grad(2), fval_grad, time_fmincon_grad);    
    fprintf('\t MSE: %.5e\n', mean(sqrt((sol_fmin_grad-sol_bfgs).^2)));    
    
    fig = figure("Name", sprintf("SQP - BFGS - Himmelblau - Solution for x0=[%+.1f, %+.1f]",x0(1),x0(2)), 'Position', [150, 150, 600, 600]);
    hold on
    % Create contour plot and constraints
    [cfig, conFigs] = contourHimmel(true);
    
    % Add points of interest
    h = traceIterations(output.xk, "b");
    intp = plotPoint(x0(1),x0(2), "int");
    solp_bfgs = plotPoint(sol_bfgs(1),sol_bfgs(2), "sol", "y", 16);
    solp_grad = plotPoint(sol_fmin_grad(1),sol_fmin_grad(2), "gen", "g", 8);
    legend([intp,solp_grad,solp_bfgs, h],{'x_0', 'fmincon solution', 'SQP BFGS sol.', "Trace for BFGS"},'Location','southwest')
    hold off
    savefigpdf(fig, sprintf("ex4_6_bfgs_himmelblau_x0=%+.0f_%+.0f",x0(1),x0(2)), 4);

    data(j,:) = [fval_grad, sol_fmin_grad', time_fmincon_grad, nan ...
                 obj, sol_bfgs', t_bfgs_total, output.iterations, ...
                 nan, mean(sqrt((sol_fmin_grad-sol_bfgs).^2))];

end

input.data = data';
% Set column labels (use empty string for no label):
input.tableColLabels = sprintfc("$x_0=[%.1f, %.1f]$", x0s')';
% Set row labels (use empty string for no label):
input.tableRowLabels = {'$f(x)_{\textit{fmincom}} =$', ...
                        '$x_{\textit{fmincom}} = $', ...
                        '', ...
                        '$\text{time}_{\textit{fmincom}}\, [s] =$', ...
                        '', ...
                        '$f(x)_{\textit{SQP BFGS}} =$', ...
                        '$x_{\textit{SQP BFGS}} =$', ...
                        '', ...
                        '$\text{time}_{\textit{SQP BFGS}}\, [s] =$', ...
                        '$\text{Interations}_{\textit{SQP BFGS}}\, =$', ...
                        '', ...
                        'MSE = '};
                        
% Set the row format of the data values 
input.dataFormatMode = 'row';
input.dataFormat = {'%.5f', 9, "%d", 1, "%.5e",2};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'r';
% Switch table borders on/off:
input.booktabs = 1;
% LaTex table caption:
input.tableCaption = sprintf('Comparison of found solution of SQP BFGS and $\textit{fmincon}$ for different inital points for the Himmelblau test problem.');
% LaTex table label:
input.tableLabel = 'ex4_bfgs_himmel';
input.makeCompleteLatexDocument = 0;
input.dataNanString = '';
input.tablePlacement = '!ht';
% Now call the function to generate LaTex code:
latex = latexTable(input);
savelatexTable(latex, input.tableLabel, 4);

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
