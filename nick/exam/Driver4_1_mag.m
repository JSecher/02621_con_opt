%% Problem 4 Driver
%% Clean up
clear
close all

runContourPlot_44 = false;
runSolveTest_45= false;
runSolveTest_452= false;
runBFGS_46 = true;
runBFGS_LS_47 = false;
himmelblau49 = true;

addpath('casadi-windows-matlabR2016a-v3.5.5')
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
%savefigpdf(fig, "ex4_himmelblau_poi", 4);

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
    %solp = plotPoint(sol_fmin(1),sol_fmin(2), "sol", 'y', 18);
    solp_grad = plotPoint(sol_fmin_grad(1),sol_fmin_grad(2), "min", 'r', 10);
    solp_cas = plotPoint(sol_cas(1),sol_cas(2), "gen", 'g', 8);
    legend([intp,solp_grad,solp_cas],{'Initial Point', 'fmincon solution','CasADi solution'})
    hold off
    %savefigpdf(fig, sprintf("ex4_5_himmelblau_x0=%+.0f_%+.0f",x0(1),x0(2)), 4);

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
%savelatexTable(latex, input.tableLabel, 4);


end

%% Problem 4.5 - Solve Rosenbrock Part 1
if runSolveTest_452
%%
disp("Solving Rosenbrock Test problem using fmincon and Casadi");

% Select constraint for rosenbrok, 1 = circle, 2 = box
con_case = 2;

switch con_case
    case 1
        con_type_name = "Circle";
        x0s = [[0.0; 0.0]];
    case 2
        con_type_name = "Box";
        x0s = [[-0.5; 1]];
    otherwise
        error("Unknown case for rosenbrock")
end

fig = figure("Name", sprintf("Rosenbrock - solutions for x0 %s constraints", con_type_name), 'Position', [100, 100, 800,800]);
hold on
clear data
% Create contour plot and constraints
[cfig, conFigs] = contourRosen(con_case);

for j=1:1
    % Define variables for casadi
    x1 = SX.sym('x1');  % Define variables
    x2 = SX.sym('x2');
    
    % Solve problem using fmincon
    options = optimoptions('fmincon',... 
                           'Display','none',... 
                           'Algorithm','interior-point');

    switch con_case
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
            
            % For casadi
            g = [x1^2+x2^2];
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
            % For casadi
            g = [];
        otherwise
            error("Unknown case for rosenbrock")
    end

    % Define problem
    hbprob = struct('x',[x1; x2], ...
                    'f',100 * (x2-x1^2)^2 + (1-x1)^2, ...
                    'g',g);
    
    % Solve the problem with casadi
    options = struct;
    options.ipopt.print_level = 0;
    options.print_time = 0;
    tstart = cputime;
    S = nlpsol('S', 'ipopt', hbprob, options);
    switch con_case
        case 1
            r = S('x0', x0, 'lbg',0,'ubg',1, 'lbx',-1,'ubx',1);
        case 2
            r = S('x0', x0, 'lbx',[-1.5, 0.5],'ubx',[1.5, 1.5]);
    end
    time_cas = cputime - tstart;
    sol_cas = full(r.x);
    obj_cas = full(r.f);
    
    % Print results
    fprintf('Found solutions for x0 = [%.2f, %.2f] ::\n', x0(1), x0(2)); 
    fprintf('\t fmincon solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_fmin(1), sol_fmin(2), fval, time_fmincon);
    fprintf('\t Casadi solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_cas(1), sol_cas(2), obj_cas, time_cas);
    
    data(j,:) = [fval, sol_fmin', time_fmincon, nan, ...
                 obj_cas, sol_cas', time_cas, nan, mean(sqrt((sol_fmin-sol_cas).^2))];

        
    % Add points of interest
    %fm = traceIterations([x0, sol_fmin], "r", '-');
    %cas = traceIterations([x0, sol_cas], "g");
    intp = plotPoint(x0(1),x0(2), "int");
    solp_fm = plotPoint(sol_fmin(1),sol_fmin(2), "sol", "r", 18);
    solp_cas = plotPoint(sol_cas(1),sol_cas(2), "gen", "g", 8);
end

legend([intp,solp_fm,solp_cas, conFigs(1)],{'Initial Point', 'fmincon solution', 'CasADi solution',})
hold off
%savefigpdf(fig, sprintf('ex4_solve_x0s_test_rosen_%s', lower(con_type_name)), 4);

input.data = data';
% Set column labels (use empty string for no label):
input.tableColLabels = sprintfc("$x_0=[%.1f, %.1f]$", x0s')';
% Set row labels (use empty string for no label):
input.tableRowLabels = {'$f(x)_{\textit{fmincom}} =$', ...
                        '$x_{\textit{fmincom}} = $', '', ...
                        '$\text{time}_{\textit{fmincom}}\, [s] =$', ...
                        '', ...
                        '$f(x)_{\textit{CasADi}} =$', ...
                        '$x_{\textit{CasADi}} = $', '', ...
                        '$\text{time}_{\textit{CasADi}}\, [s] =$', ...
                        '', 'Mean Squared Difference = '};
% Set the row format of the data values 
%input.dataFormatMode = 'row';
input.dataFormat = {'%.4f'};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'r';
% Switch table borders on/off:
input.booktabs = 1;
% LaTex table caption:
input.tableCaption = sprintf('Found solution for different inital points for the Rosenbrock test problem with %s constraint .', lower(con_type_name));
% LaTex table label:
input.tableLabel = sprintf('ex4_solve_test_rosen_%s', lower(con_type_name));
input.makeCompleteLatexDocument = 0;
input.dataNanString = '';
input.tablePlacement = '!ht';
% Now call the function to generate LaTex code:
latex = latexTable(input);
%savelatexTable(latex, input.tableLabel, 4);

%% Contour plot - Rosenbrock
f = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;
view = 6;
res = 50;
x = -view:1/res:view;
y = -view:1/res:view;
[X,Y] = meshgrid(x,y);
v = [0:2:10 10:10:100 100:20:300];
F = f(X,Y);
fig = figure;
hold on
[m,c]=contour(X,Y,F,v,"linewidth",2);
colormap("turbo")
xlabel('x_1','Fontsize',14)
ylabel('x_2','Fontsize',14)
colorbar
xlim([-2, 2])
ylim([-2, 4])
hold off
%savefigpdf(fig, "ex4_5_rosenbrock", 4);
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
    %solp_grad = plotPoint(sol_fmin_grad(1),sol_fmin_grad(2), "gen", "g", 8);
    legend([intp,solp_bfgs, h],{'Initial Point', 'SQP BFGS sol.', "Trace for BFGS"},'Location','southwest')
    hold off
    %savefigpdf(fig, sprintf("ex4_6_bfgs_himmelblau_x0=%+.0f_%+.0f",x0(1),x0(2)), 4);

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
input.dataFormat = {'%.4f', 9, "%d", 1, "%.4e",2};
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
%savelatexTable(latex, input.tableLabel, 4);

end
%% Problem 4.7 - Dampend BFGS with line search 
if runBFGS_LS_47

clear data
x0s = [[0.0; 0.0], [-4.0; 0], [-4; 1]];

fig = figure("Name", "SQP - Line Search - Himmelblau - Solution for x0s", 'Position', [150, 150, 600, 600]);
hold on
% Create contour plot and constraints
[cfig, conFigs] = contourHimmel(true);

traces = [];
legens = [];
for j=1:length(x0s)
    %sympref('FloatingPointOutput',1);
    x0 = x0s(:,j);    % Initial point
    
    xl = [-5; -5];      % Lower bound for x
    xu = [5; 5];        % Upper bound for x
    cl = [0; 0];        % Lower bound for constraints 
    cu = [47; 70];      % Upper bound for constraints
    
    tstart = cputime;
    [sol_ls, obj, lambda, output] = SQPsolver(@objfungradHimmelblau, ...
                                                @confungradHimmelblau1, ...
                                                xl, xu, ...
                                                cl, cu, ...
                                                x0, 'line');
    t_ls_total = cputime - tstart;
    
    
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
    fprintf('\t Line Search solution: [%.5e, %.5e], objective: %.5e, time: %.5e, iter: %d\n', sol_ls(1), sol_ls(2), obj, t_ls_total, output.iterations);
    fprintf('\t fmincon grad solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_fmin_grad(1), sol_fmin_grad(2), fval_grad, time_fmincon_grad);    
    fprintf('\t MSE: %.5e\n', mean(sqrt((sol_fmin_grad-sol_ls).^2)));    
    

    % Add points of interest
    h = traceIterations(output.xk, "b");
    intp = plotPoint(x0(1),x0(2), "int");
    solp_ls = plotPoint(sol_ls(1),sol_ls(2), "sol", "y", 16);
    % solp_grad = plotPoint(sol_fmin_grad(1),sol_fmin_grad(2), "gen", "g", 8);
    
    data(j,:) = [fval_grad, sol_fmin_grad', time_fmincon_grad, nan ...
                obj, sol_ls', t_ls_total, output.iterations, output.function_calls ...
                 nan, mean(sqrt((sol_fmin_grad-sol_ls).^2))];

end

legend([intp,solp_ls, h],{'Initial Point' 'SQP LS sol.', 'Trace for LS.'},'Location','southwest')
hold off
%savefigpdf(fig, "ex4_6_ls_himmelblau_x0s", 4);


input.data = data';
% Set column labels (use empty string for no label):
input.tableColLabels = sprintfc("$x_0=[%.1f, %.1f]$", x0s')';
% Set row labels (use empty string for no label):
input.tableRowLabels = {'$f(x)_{\textit{fmincom}} =$', ...
                        '$x_{\textit{fmincom}} = $', ...
                        '', ...
                        '$\text{time}_{\textit{fmincom}}\, [s] =$', ...
                        '', ...
                        '$f(x)_{\textit{SQP LS}} =$', ...
                        '$x_{\textit{SQP LS}} =$', ...
                        '', ...
                        '$\text{time}_{\textit{SQP LS}}\, [s] =$', ...
                        '$\text{Interations}_{\textit{SQP LS}}\, =$', ...
                        '$\text{Function Calls}_{\textit{SQP LS}}\, =$',...
                        '', ...
                        'MSE = '};
                        
% Set the row format of the data values 
input.dataFormatMode = 'row';
input.dataFormat = {'%.5f', 9, "%d", 2, "%.5e",2};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'r';
% Switch table borders on/off:
input.booktabs = 1;
% LaTex table caption:
input.tableCaption = sprintf('Comparison of found solution of SQP Line Search and \\textit{fmincon} for different inital points for the Himmelblau test problem.');
% LaTex table label:
input.tableLabel = 'ex4_ls_himmel';
input.makeCompleteLatexDocument = 0;
input.dataNanString = '';
input.tablePlacement = '!ht';
% Now call the function to generate LaTex code:
latex = latexTable(input);
%savelatexTable(latex, input.tableLabel, 4);

end


%% 4.9 Himmelblau
if himmelblau49
clear data
x0s = [[-4.0; 0]];

fig = figure("Name", "SQP - Line Search - Himmelblau - Solution for x0s", 'Position', [150, 150, 600, 600]);
hold on
% Create contour plot and constraints
[cfig, conFigs] = contourHimmel(true);

traces = [];
legens = [];
for j=1:1
    %sympref('FloatingPointOutput',1);
    x0 = x0s(:,j);    % Initial point
    
    xl = [-5; -5];      % Lower bound for x
    xu = [5; 5];        % Upper bound for x
    cl = [0; 0];        % Lower bound for constraints 
    cu = [47; 70];      % Upper bound for constraints
    
    tstart = cputime;
    [sol_ls, obj, lambda, output] = SQPsolver(@objfungradHimmelblau, ...
                                                @confungradHimmelblau1, ...
                                                xl, xu, ...
                                                cl, cu, ...
                                                x0, 'line');
    t_ls_total = cputime - tstart;
    
    
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

    fprintf('Found solutions for x0 = [%.2f, %.2f] ::\n', x0(1), x0(2)); 
    fprintf('\t Line Search solution: [%.5e, %.5e], objective: %.5e, time: %.5e, iter: %d\n', sol_ls(1), sol_ls(2), obj, t_ls_total, output.iterations);
    fprintf('\t fmincon grad solution: [%.5e, %.5e], objective: %.5e, time: %.5e\n', sol_fmin_grad(1), sol_fmin_grad(2), fval_grad, time_fmincon_grad);    
    fprintf('\t MSE: %.5e\n', mean(sqrt((sol_fmin_grad-sol_ls).^2)));    
    

    % Add points of interest
    h = traceIterations(output.xk, "b");
    intp = plotPoint(x0(1),x0(2), "int");
    solp_ls = plotPoint(sol_ls(1),sol_ls(2), "sol", "y", 16);
    solp_grad = plotPoint(sol_fmin_grad(1),sol_fmin_grad(2), "gen", "g", 8);
    solp_cas = plotPoint(sol_cas(1),sol_cas(2), "min", "r", 10);
    
    data(j,:) = [fval_grad, sol_fmin_grad', time_fmincon_grad, nan ...
                obj, sol_ls', t_ls_total, output.iterations, output.function_calls ...
                 nan, mean(sqrt((sol_fmin_grad-sol_ls).^2))];

end

legend([intp, solp_grad, solp_cas, solp_ls, h],{'Initial Point', 'fmincon sol.', 'CasADI sol', 'SQP LS sol.', 'Trace for LS.'},'Location','southwest')
hold off
%savefigpdf(fig, "ex4_6_ls_himmelblau_x0s", 4);

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

function f = objfunRosenbrock(x,p)
    % Function giving the objective of the Rosenbrock problem
    f = 100 * (x(2)-x(1)^2)^2 + (1-x(1))^2;
end

function [c,ceq] = confunRosenbrock(x,p)
    % Function giving the constraints for the Himmelblau problem
    % Function taken from: Slide 8, Lecture 01B 
    c = zeros(1,1);
    ceq = zeros(0,1);

    % Inequality constraints c(x) <= 0
    c(1,1) = x(1)^2 + x(2)^2 - 1;
end

function [c,ceq] = confunRosenbrock_part2(x,p)
    % Function giving the constraints for the Himmelblau problem
    % Function taken from: Slide 8, Lecture 01B 
    c = zeros(0,1);
    ceq = zeros(0,1);
end

function [f,dfdx] = objfungradRosenbrock(x,p)
    % Function giving the objective and gradients of the Himmelblau problem
    % Function taken from: Slide 8, Lecture 01B 
    f = 100 * (x(2)-x(1)^2)^2 + (1-x(1))^2;
    % compute the gradient of f
    if nargout > 1
        dfdx = zeros(2,1);
        dfdx(1,1) = -400*(-x(1)^2 + x(2))*x(1) - 2 + 2*x(1);
        dfdx(2,1) = -200*x(1)^2 + 200*x(2);
    end
end

function [c,ceq,dcdx,dceqdx] = confungradRosenbrock(x,p) 
    % Function giving the constraints and gradients of the Himmelblau problem
    % Function taken from: Slide 8, Lecture 01B 
    c = zeros(1,1);
    ceq = zeros(0,1);
    % Inequality constraints c(x) <= 0 
    c(1,1) = x(1)^2 + x(2)^2 - 1;
    % Compute constraint gradients
    if nargout > 2
        dcdx = zeros(2,1); 
        dceqdx = zeros(2,0);
        dcdx(1,1) = 2*x(1); % dc1dx1
        dcdx(2,1) = 2*x(2); % dc1dx2
    end
end

function [c,ceq,dcdx,dceqdx] = confungradRosenbrock_part2(x,p) 
    % Function giving the constraints and gradients of the Himmelblau problem
    % Function taken from: Slide 8, Lecture 01B 
    c = zeros(0,1);
    ceq = zeros(0,1);
    % Inequality constraints c(x) <= 0 
    % Compute constraint gradients
    if nargout > 2
        dcdx = zeros(2,1); 
        dceqdx = zeros(2,0);
    end
end