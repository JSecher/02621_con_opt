%% Problem 3.4 - 3.6 Driver
%% Clean up
clear
close all

% Options for linprog
default_options = optimset('Display', 'off');

run_test1 = false;
run_test2 = true;
run_test3 = true;

%% Test on toy problem for correctness
if run_test1

Nsamples = 10; % Run test multiple times for small problems
% Defining the data
g = [-16.1000,  -8.5000,  -15.7000, -10.0200, -18.6800]';

A = ones(5,1);
b = ones(1,1);

[n_var,m_eqcon] = size(A);

l = zeros(n_var,1);
u = ones(n_var,1);

% Starting point 
x0 = zeros(n_var,1);
s0 = ones(2*n_var,1);
y0 = ones(m_eqcon,1);
z0 = ones(2*n_var,1);


% Test on own, Linprog interior point and Linprog simplex

% Test on own
tic;
for i=1:Nsamples
[x,y,z,s,iter] = LPSolverInteriorPoint(g,A,b,l,u,x0,y0,z0,s0);
end
time_own = toc/Nsamples;
obj = g'*x;

% Test Linprog interior-point
options = optimset(default_options, 'Algorithm', 'interior-point');
tic;
for i=1:Nsamples
[x_linprog, obj_linprog, ~, output] = linprog(g,[],[],A',b,l,u, options);
end
time_linprog = toc/Nsamples;
iter_linprog = output.iterations;

% Test Linprog Dual Simplex
options = optimset(default_options, 'Algorithm', 'dual-simplex');
tic;
for i=1:Nsamples
[x_simplex, obj_simplex, ~, output] = linprog(g,[],[], A', b,l,u, options);
end
time_simplex= toc/Nsamples;
iter_simplex = output.iterations;

err = mean(sqrt((x - x_simplex).^2));
err_linprog = mean(sqrt((x_simplex - x_linprog).^2));

fprintf("Linprog Algo | Obj: %.5f \t Time: %.5f \t Iter: %3d \t MSE %.5e, obj Err %5e \n", obj_linprog, time_linprog, iter_linprog, err_linprog, obj_linprog-obj_simplex)
fprintf("Own Algo     | Obj: %.5f \t Time: %.5f \t Iter: %3d \t MSE %.5e, obj Err %5e \n", obj, time_own, iter, err, obj-obj_simplex)
fprintf("Simplex      | Obj: %.5f \t Time: %.5f \t Iter: %3d\n", obj_simplex, time_simplex, iter_simplex)

fprintf("\n Thats nice...!, lets make some latex out of that, \n \n * *abrakadabra*\n \n *pafff*\n")

%% Make Table


data = [obj_linprog, obj];
data(2,:) = [nan, obj-obj_linprog];
data(3:7,:) = [x_linprog, x];
data(8,:) = [nan, err];
data(9,:) = [nan, nan];
data(10,:) = [time_linprog, time_own];
data(11,:) = [iter_linprog, iter];
input.data = data;

% Set column labels (use empty string for no label):
input.tableColLabels = {'Linprog IP', 'Own Solver'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'$\phi =$', ...
                        '$\Delta \phi =$', ...
                        'x =', '', '', '', '', ...
                        'MSE =', ...
                        '', ...
                        'Iterations =','Time ='};
% Set the row format of the data values 
% Set the row format of the data values 
input.dataFormatMode = 'row';
input.dataFormat = {'%.4e',10,'%d',1};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'r';
% Switch table borders on/off:
input.booktabs = 1;
% LaTex table caption:
input.tableCaption = sprintf('Comparison between LinProg interior-point and Own LP interior-point alogorithm. Error is computed relative to LinProg IP solver.');
% LaTex table label:
input.tableLabel = 'ex3_objective_comp';
input.makeCompleteLatexDocument = 0;
input.dataNanString = '';
input.tablePlacement = 'ht';
% Now call the function to generate LaTex code:
latex = latexTable(input);
savelatexTable(latex, input.tableLabel, 3);

% Make table with simplex also

data = [obj_linprog, obj, obj_simplex];
data(2,:) = [obj_linprog-obj_simplex, obj-obj_simplex, nan];
data(3:7,:) = [x_linprog, x, x_simplex];
data(8,:) = [err_linprog, err, nan];
data(9,:) = [nan, nan, nan];
data(10,:) = [time_linprog, time_own, time_simplex];
data(11,:) = [iter_linprog, iter, iter_simplex];
input.data = data;


% Set column labels (use empty string for no label):
input.tableColLabels = {'Linprog IP', 'Own Solver', 'Linprog Simplex'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'$\phi =$', ...
                        '$\Delta \phi =$', ...
                        'x =', '', '', '', '', ...
                        'MSE =', ...
                        '', ...
                        'Iterations =','Time ='};
% Set the row format of the data values 
input.dataFormatMode = 'row';
input.dataFormat = {'%.4e',10,'%d',1};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'r';
% Switch table borders on/off:
input.booktabs = 1;
% LaTex table caption:
input.tableCaption = sprintf('Comparison between LinProg interior-point, Own LP interior-point and LinProg Simplex alogorithm. Error is computed relative to LinProg simplex solver.');
% LaTex table label:
input.tableLabel = 'ex3_objective_comp_simplex';
input.makeCompleteLatexDocument = 0;
input.dataNanString = '';
input.tablePlacement = 'ht';
% Now call the function to generate LaTex code:
latex = latexTable(input);
savelatexTable(latex, input.tableLabel, 3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Size dependent problem 1    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_test2
% Number of tests and test sizes and arrays for storing results
Ntests = 20;

Nsizes = round(linspace(100,2000, Ntests));

solvers = ["Own Solver", "LinProg IP", "LinProg Simplex"];
Nsolvers = length(solvers);

objs = zeros(Nsolvers,Ntests);
iterations = zeros(Nsolvers,Ntests);
times = zeros(Nsolvers, Ntests);
MSE_errs = zeros(Nsolvers-1,Ntests);
obj_errs = zeros(Nsolvers-1,Ntests);

x_true = nan*zeros(Ntests,max(Nsizes));

fprintf("\nTest 2 : Size dendendt problem, %d solvers, for each %d value(s)" + ...
        " between [%d, %d].\n\n", Nsolvers, Ntests,Nsizes(1),Nsizes(end));

rng(200) % Set seed
poorMansProgressBar(Ntests);
for j=1:Ntests
    poorMansProgressBar(0);  % Update progress for iteration
        
    
    n_val = Nsizes(j);  % Set value of b
    
    % Generate problem
    % Defining the data
    [g,A,b,x,lam] = randomLP(n_val,floor(n_val/2), 1);
    [n_var,m_eqcon] = size(A);
    l = zeros(n_val,1);
    u = ones(n_val,1);

    % Starting point 
    x0 = zeros(n_var,1);
    s0 = ones(2*n_var,1);
    y0 = ones(m_eqcon,1);
    z0 = ones(2*n_var,1);

    % Run Own 
    tic;
    [x,y,z,s,iter] = LPSolverInteriorPoint(g,A,b,l,u,x0,y0,z0,s0);
    time = toc;
    obj = g'*x;

    % Run Linprog IP 
    options = optimset(default_options, 'Algorithm', 'interior-point');
    tic;
    [x_ip, obj_ip, ~, outtmp] = linprog(g,[],[],A',b,l,u, options);
    time_ip = toc;
    iter_ip = outtmp.iterations;
    
    % Run Linprog Simplex
    options = optimset(default_options, 'Algorithm', 'dual-simplex');
    tic;
    [x_sim, obj_sim, ~, outtmp] = linprog(g,[],[],A',b,l,u, options);
    time_sim = toc;
    iter_sim = outtmp.iterations;
    
    % Save and compute stuff
    objs(1,j) = obj;
    objs(2,j) = obj_ip;
    objs(3,j) = obj_sim;

    obj_errs(1,j) = obj - obj_ip;
    obj_errs(2,j) = obj - obj_sim;

    MSE_errs(1,j) = mean(sqrt((x - x_ip).^2));
    MSE_errs(2,j) = mean(sqrt((x - x_sim).^2));
    
    iterations(1,j) = iter;
    iterations(2,j) = iter_ip;
    iterations(3,j) = iter_sim;
    times(1,j) = time;
    times(2,j) = time_ip;
    times(3,j) = time_sim;
    
end
poorMansProgressBar(-1);

%% Make Plots for large problem


f = figure('Name','Objective');
for i=1:Nsolvers
    plot(Nsizes, objs(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of variables, n")
ylabel("Objective value")
savefigpdf(f, "ex3_obj_large", 3);

f = figure('Name','Objective semilogy');
for i=1:Nsolvers
    semilogy(Nsizes, objs(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of variables, n")
ylabel("Objective value")
savefigpdf(f, "ex3_obj_semilogy_large", 3);

f = figure('Name','Time');
for i=1:Nsolvers
    plot(Nsizes, times(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of variables, n")
ylabel("Time [s]")
savefigpdf(f, "ex3_time_large", 3);

f = figure('Name','Time log');
for i=1:Nsolvers
    semilogy(Nsizes, times(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of variables, n")
ylabel("Time [s]")
savefigpdf(f, "ex3_time_semilogy_large", 3);

f = figure('Name','Iterations');
for i=1:Nsolvers
    plot(Nsizes, iterations(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of variables, n")
ylabel("Iterations")
savefigpdf(f, "ex3_iterations_large", 3);

f = figure('Name','Iterations semilogy');
for i=1:Nsolvers
    semilogy(Nsizes, iterations(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of variables, n")
ylabel("Iterations")
savefigpdf(f, "ex3_iterations_semilogy_large", 3);


f = figure('Name','Iterations No simplex');
for i=1:Nsolvers-1
    semilogy(Nsizes, iterations(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers([1,2]), 'Location','northwest')
xlabel("Number of variables, n")
ylabel("Iterations")
savefigpdf(f, "ex3_iterations_no_simplex_large", 3);

f = figure('Name','Objective Error');
for i=1:Nsolvers-1
    plot(Nsizes, obj_errs(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(sprintfc("Obj. err. vs %s",solvers(2:end)), 'Location','southeast')
xlabel("Number of variables, n")
ylabel("Objective error")
savefigpdf(f, "ex3_obj_err_large", 3);

f = figure('Name','Abs Objective Error semilog');
for i=1:Nsolvers-1
    semilogy(Nsizes, abs(obj_errs(i,:)), 'LineWidth',2);
    hold on
end
grid on
legend(sprintfc("Obj. err. vs %s",solvers(2:end)), 'Location','northwest')
xlabel("Number of variables, n")
ylabel("Absolute Objective error")
savefigpdf(f, "ex3_obj_err_semilogy_large", 3);

f = figure('Name','MSE Error');
for i=1:Nsolvers-1
    plot(Nsizes, MSE_errs(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(sprintfc("MSE vs %s",solvers(2:end)), 'Location','northwest')
xlabel("Number of variables, n")
ylabel("MSE")
savefigpdf(f, "ex3_mse_err_large", 3);

f = figure('Name','MSE Error log');
for i=1:Nsolvers-1
    semilogy(Nsizes, MSE_errs(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(sprintfc("MSE vs %s",solvers(2:end)), 'Location','northwest')
xlabel("Number of variables, n")
ylabel("MSE")
savefigpdf(f, "ex3_mse_err_semilogy_large", 3);


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Size dependent problem 2 - Growing constraints    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_test3
% Number of tests and test sizes and arrays for storing results
Ntests = 20;

Nfix = 2000;
Msizes = round(linspace(100,Nfix, Ntests), -1);

solvers = ["Own Solver", "LinProg IP", "LinProg Simplex"];
Nsolvers = length(solvers);

objs = zeros(Nsolvers,Ntests);
iterations = zeros(Nsolvers,Ntests);
times = zeros(Nsolvers, Ntests);
MSE_errs = zeros(Nsolvers-1,Ntests);
obj_errs = zeros(Nsolvers-1,Ntests);


fprintf("\nTest 3 : Size dendendt problem - growing constraints, n fixed to n=%d, %d solvers, for each %d number of constraint(s)" + ...
        " between [%d, %d].\n\n", Nfix, Nsolvers, Ntests,Msizes(1),Msizes(end));

rng(200) % Set seed
poorMansProgressBar(Ntests);
for j=1:Ntests
    poorMansProgressBar(0);  % Update progress for iteration
        
    
    m_val = Msizes(j);  % Set value of b
    
    % Generate problem
    % Defining the data
    [g,A,b,x,lam] = randomLP(Nfix,m_val, 1);
    [n_var,m_eqcon] = size(A);
    l = zeros(Nfix,1);
    u = ones(Nfix,1);

    % Starting point 
    x0 = zeros(Nfix,1);
    s0 = ones(2*Nfix,1);
    y0 = ones(m_val,1);
    z0 = ones(2*Nfix,1);

    % Run Own 
    tic;
    [x,y,z,s,iter] = LPSolverInteriorPoint(g,A,b,l,u,x0,y0,z0,s0);
    time = toc;
    obj = g'*x;

    % Run Linprog IP 
    options = optimset(default_options, 'Algorithm', 'interior-point');
    tic;
    [x_ip, obj_ip, ~, outtmp] = linprog(g,[],[],A',b,l,u, options);
    time_ip = toc;
    iter_ip = outtmp.iterations;
    
    % Run Linprog Simplex
    options = optimset(default_options, 'Algorithm', 'dual-simplex');
    tic;
    [x_sim, obj_sim, ~, outtmp] = linprog(g,[],[],A',b,l,u, options);
    time_sim = toc;
    iter_sim = outtmp.iterations;
    
    % Save and compute stuff
    objs(1,j) = obj;
    objs(2,j) = obj_ip;
    objs(3,j) = obj_sim;

    obj_errs(1,j) = obj - obj_ip;
    obj_errs(2,j) = obj - obj_sim;

    MSE_errs(1,j) = mean(sqrt((x - x_ip).^2));
    MSE_errs(2,j) = mean(sqrt((x - x_sim).^2));
    
    iterations(1,j) = iter;
    iterations(2,j) = iter_ip;
    iterations(3,j) = iter_sim;
    times(1,j) = time;
    times(2,j) = time_ip;
    times(3,j) = time_sim;
    
end
poorMansProgressBar(-1);

%% Make Plots for large problem


f = figure('Name','Growing m - Objective');
for i=1:Nsolvers
    plot(Msizes, objs(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("Objective value")
savefigpdf(f, "ex3_obj_large_constraints", 3);

f = figure('Name','Growing m - Time');
for i=1:Nsolvers
    plot(Msizes, times(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("Time [s]")
savefigpdf(f, "ex3_time_large_constraints", 3);

f = figure('Name','Growing m - Time log');
for i=1:Nsolvers
    semilogy(Msizes, times(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("Time [s]")
savefigpdf(f, "ex3_time_semilogy_large_constraints", 3);

f = figure('Name','Growing m - Iterations');
for i=1:Nsolvers
    plot(Msizes, iterations(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("Iterations")
savefigpdf(f, "ex3_iterations_large_constraints", 3);

f = figure('Name','Growing m - Iterations semilog');
for i=1:Nsolvers
    semilogy(Msizes, iterations(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers, 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("Iterations")
savefigpdf(f, "ex3_iterations_semilog_large_constraints", 3);


f = figure('Name','Growing m - Iterations No simplex');
for i=1:Nsolvers-1
    plot(Msizes, iterations(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(solvers([1,2]), 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("Iterations")
savefigpdf(f, "ex3_iterations_no_sim_large_constraints", 3);

f = figure('Name','Growing m - Objective Error');
for i=1:Nsolvers-1
    plot(Msizes, obj_errs(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(sprintfc("Obj. err. vs %s",solvers(2:end)), 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("Objective error")
savefigpdf(f, "ex3_obj_err_large_constraints", 3);

f = figure('Name','Growing m - abs Objective Error semilog');
for i=1:Nsolvers-1
    semilogy(Msizes, abs(obj_errs(i,:)), 'LineWidth',2);
    hold on
end
grid on
legend(sprintfc("Obj. err. vs %s",solvers(2:end)), 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("Absolute Objective error")
savefigpdf(f, "ex3_obj_err_semilogy_large_constraints", 3);

f = figure('Name','Growing m - MSE Error Semilog');
for i=1:Nsolvers-1
    plot(Msizes, MSE_errs(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(sprintfc("MSE vs %s",solvers(2:end)), 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("MSE")
savefigpdf(f, "ex3_mse_err_large_constraints", 3);

f = figure('Name','Growing m - MSE Error Semilog');
for i=1:Nsolvers-1
    semilogy(Msizes, MSE_errs(i,:), 'LineWidth',2);
    hold on
end
grid on
legend(sprintfc("MSE vs %s",solvers(2:end)), 'Location','northwest')
xlabel("Number of constraints, m")
ylabel("MSE")
savefigpdf(f, "ex3_mse_err_semilogy_large_constraints", 3);


end


%% Nothing good below here

function [msg] = poorMansProgressBar(state)
% state > 0 : Begin progress bar, state gives expected number of iterations
% state = 0 : Update progress bar, call this max the expected number of
% times
% state = -1 : Display done
if state > 0
    msg = sprintf("Start |");
    for i=1:state                 
        msg = msg + sprintf(" -");
    end                             
    msg = msg + sprintf(" | Finish\n       ");
elseif state == 0
    msg = sprintf(" *");
else
    msg = fprintf("   DONE!!\n");
end
fprintf("%s",msg);
end



