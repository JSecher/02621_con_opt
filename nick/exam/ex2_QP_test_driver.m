%% Problem 2.4 - 2.6 Driver
%% Clean up
clear
close all

%%
% Defining the data
H = [ 5.0000,    1.8600,    1.2400,   1.4800   -0.4600;
      1.8600,    3.0000,    0.4400,   1.1200,   0.5200;
      1.2400,    0.4400,    3.8000,   1.5600,  -0.5400;
      1.4800,    1.1200,    1.5600,   7.2000,  -1.1200;
     -0.4600,    0.5200,   -0.5400,  -1.1200,   7.8000];

g = [-16.1000,  -8.5000,  -15.7000, -10.0200, -18.6800]';

A = [16.1000, 1.0000;
      8.5000, 1.0000;
     15.7000, 1.0000;
     10.0200, 1.0000;
     18.6800, 1.0000];

b = [15, 1]';

[n_var,m_eqcon] = size(A);
l = zeros(n_var,1);
u = ones(n_var,1);

% Starting point 
x0 = zeros(n_var,1);
s0 = ones(2*n_var,1);
y0 = ones(m_eqcon,1);
z0 = ones(2*n_var,1);

% Generate values for b(1) in the range [8.5 18.68]
Ntests = 100;
Nsamples = 5;  % For small tests, run the function multiple times and take mean
bs = linspace(8.5, 18.68, Ntests);

% Simple test vs quadprog
xs = zeros(Ntests, n_var);
xs_qp = zeros(Ntests, n_var);
objs = zeros(Ntests,1);
objs_qp = zeros(Ntests,1);
iterations = zeros(Ntests,1);
iterations_qp = zeros(Ntests,1);
times = zeros(Ntests,1);
times_qp = zeros(Ntests,1);
err = zeros(Ntests,1);


id_warning = [];
msg_warning = {};
fprintf("Test 1 : Running for %d value(s)" + ...
        " of b(1), each test run %d time(s).\n\n", Ntests, Nsamples);

poorMansProgressBar(Ntests);
for i=1:Ntests
    poorMansProgressBar(0);  % Update progress for iteration
    b(1) = bs(i);  % Set value of b
    
    % Run Quadprog
    tic
    for sample=1:Nsamples
        [xqp, objqp, qpflag, output] = quadprog(H,g,[],[],A',b, l, u, 0, optimset('Display', 'off')); 
    end
    xs_qp(i, :) = xqp; 
    times_qp(i) = toc/Nsamples;
    objs_qp(i) = objqp;
    iterations_qp(i) = output.iterations;

    % Run Own solver
    warning('off','all')
    lastwarn("","")
    tic
    for sample=1:Nsamples
        [x,y,z,s,iter] = QPSolverMehrotraInteriorPoint(H,g,A,b,l,u,x0,y0,z0,s0);
    end
    xs(i, :) = x; 
    times(i) = toc/Nsamples;
    objs(i) = 0.5*x'*H*x + g'*x;
    iterations(i) = iter;
    w = warning('query','last');
    
    % now if a warning was raised, warnMsg and warnId will not be empty.
    [warnMsg, warnId] = lastwarn();
    % you can check the warning message or id, or just throw the warning as an error if desired
    if(~isempty(warnId))
        id_warning = [id_warning, i];
        msg_warning = [msg_warning, warnMsg];
    end
    warning('on','all')

    err(i) = mean(sqrt((x - xqp).^2));

end
poorMansProgressBar(-1);

% Print the warnings if any
if ~isempty(id_warning)
    fprintf("\n!!!! %d values raise a warning for own solver !!!!!\n", length(id_warning))
    for i=1:length(id_warning)
        fprintf("Iteration %d, b(1)=%.5f, Got a warning!\n", id_warning(i), bs(id_warning(i)))
        fprintf("\t Message: %s\n", msg_warning{i})
    end
end


%% Make Plots

f = figure('Name','b vals Relative Error solver');
%f.Position(3:4) = [800, 400];
plot(bs, err, 'LineWidth',2)
grid on
legend("Own Solver")
xlabel("b(1)")
ylabel("MSE compared to quadprog")
savefigpdf(f, "ex2_error_comp_basic", 2);

f = figure('Name','b vals Relative Error solver, semilogy');
%f.Position(3:4) = [800, 400];
semilogy(bs, err, 'LineWidth',2)
grid on
legend("Own Solver")
xlabel("b(1)")
ylabel("MSE compared to quadprog")
savefigpdf(f, "ex2_error_comp_basic_semilogy", 2);

f = figure('Name','b vals Objective Error comparison');
%f.Position(3:4) = [800, 400];
plot(bs, objs-objs_qp, 'LineWidth',2)
grid on
legend("Own Solver")
xlabel("b(1)")
ylabel("Difference in obejctive vs quadprog")
savefigpdf(f, "ex2_error_objective_comp_basic", 2);

f = figure('Name','b vals Time Comparion');
% f.Position(3:4) = [800, 400];
hold on
grid on
plot(bs, times, 'LineWidth',2)
plot(bs, times_qp, 'LineWidth',2)
legend("Own Solver", "QuadProg",'Location','northwest')
xlabel("b(1)")
ylabel("time [s]")
ylim([0, max(max(times), max(times_qp))*1.1])
if Nsamples > 1
    title(sprintf("Mean over %d samples", Nsamples))
end
savefigpdf(f, "ex2_time_comp_basic", 2);


f = figure('Name','b vals Iterations Comparion');
% f.Position(3:4) = [800, 400];
hold on
grid on
plot(bs, iterations, 'LineWidth',2)
plot(bs, iterations_qp, 'LineWidth',2)
leg = ["Own Solver", "QuadProg"];
legend(leg,'Location','northwest')
xlabel("b(1)")
ylabel("Iterations")
ylim([0, max(max(iterations), max(iterations_qp))*1.1])
savefigpdf(f, "ex2_iter_comp_basic", 2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Size dependendt problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of tests and test sizes and arrays for storing results
Ntests = 20;
Nsamples = 1;  % For small tests, run the function multiple times and take mean

Nsizes = round(linspace(100,2000, Ntests), -1);

objs = zeros(Ntests,1);
objs_qp = zeros(Ntests,1);
iterations = zeros(Ntests,1);
iterations_qp = zeros(Ntests,1);
times = zeros(Ntests,1);
times_qp = zeros(Ntests,1);
err = zeros(Ntests,1);

rng(100)
id_warning = [];
msg_warning = {};
fprintf("Test 2 : Size dendendt problem %d value(s)" + ...
        " between [%d, %d], each test run %d time(s).\n\n", Ntests,Nsizes(1),Nsizes(end), Nsamples);

poorMansProgressBar(Ntests);
for i=1:Ntests
    poorMansProgressBar(0);  % Update progress for iteration
    n_val = Nsizes(i);  % Set value of b
    
    % Generate problem
    % Defining the data
    [H,g,A,b,x,lam] = randomEQP(n_val,floor(n_val/2));
    [n_var,m_eqcon] = size(A);
    l = zeros(n_var,1);
    u = ones(n_var,1);

    % Starting point 
    x0 = zeros(n_var,1);
    s0 = ones(2*n_var,1);
    y0 = ones(m_eqcon,1);
    z0 = ones(2*n_var,1);

    % Run Quadprog
    tic
    for sample=1:Nsamples
        [xqp, objqp, qpflag, output] = quadprog(H,g,[],[],A',b, l, u, 0, optimset('Display', 'off')); 
    end
    times_qp(i) = toc/Nsamples;
    objs_qp(i) = objqp;
    iterations_qp(i) = output.iterations;

    % Run Own solver
    warning('off','all')
    lastwarn("","")
    tic
    for sample=1:Nsamples
        [x,y,z,s,iter] = QPSolverMehrotraInteriorPoint(H,g,A,b,l,u,x0,y0,z0,s0);
    end
    times(i) = toc/Nsamples;
    objs(i) = 0.5*x'*H*x + g'*x;
    iterations(i) = iter;
    w = warning('query','last');
    
    % now if a warning was raised, warnMsg and warnId will not be empty.
    [warnMsg, warnId] = lastwarn();
    % you can check the warning message or id, or just throw the warning as an error if desired
    if(~isempty(warnId))
        id_warning = [id_warning, i];
        msg_warning = [msg_warning, warnMsg];
    end
    warning('on','all')

    err(i) = mean(sqrt((x - xqp).^2));

end
poorMansProgressBar(-1);

% Print the warnings if any
if ~isempty(id_warning)
    fprintf("\n!!!! %d values raise a warning for own solver !!!!!\n", length(id_warning))
    for i=1:length(id_warning)
        fprintf("Iteration %d, b(1)=%.5f, Got a warning!\n", id_warning(i), bs(id_warning(i)))
        fprintf("\t Message: %s\n", msg_warning{i})
    end
end


%% Make Plots for large problem

f = figure('Name','Large Problem b vals Relative Error solver');
f.Position(3:4) = [800, 400];
plot(Nsizes, err, 'LineWidth',2)
grid on
legend("Own Solver")
xlabel("Number of variables, n")
ylabel("MSE compared to quadprog")
savefigpdf(f, "ex2_error_comp_large", 2);

f = figure('Name','Large Problem b vals Relative Error solver, semilogy');
f.Position(3:4) = [800, 400];
semilogy(Nsizes, err, 'LineWidth',2)
grid on
legend("Own Solver")
xlabel("Number of variables, n")
ylabel("MSE compared to quadprog")
savefigpdf(f, "ex2_error_comp_large_semilogy", 2);

f = figure('Name','Large Problem b vals Objective Error comparison');
%f.Position(3:4) = [800, 400];
plot(Nsizes, objs-objs_qp, 'LineWidth',2)
grid on
legend("Own Solver")
xlabel("b(1)")
ylabel("Difference in obejctive vs quadprog")
savefigpdf(f, "ex2_error_objective_comp_large", 2);

f = figure('Name','Large Problem b vals Time Comparion');
% f.Position(3:4) = [800, 400];
hold on
grid on
plot(Nsizes, times, 'LineWidth',2)
plot(Nsizes, times_qp, 'LineWidth',2)
legend("Own Solver", "QuadProg")
xlabel("Number of variables, n")
ylabel("time [s]")
ylim([0, max(max(times), max(times_qp))*1.1])
if Nsamples > 1
    title(sprintf("Mean over %d samples", Nsamples))
end
savefigpdf(f, "ex2_time_comp_large", 2);


f = figure('Name','Large Problem b vals Iterations Comparion');
% f.Position(3:4) = [800, 400];
hold on
grid on
plot(Nsizes, iterations, 'LineWidth',2)
plot(Nsizes, iterations_qp, 'LineWidth',2)
legend("Own Solver", "QuadProg")
xlabel("Number of variables, n")
ylabel("Iterations")
ylim([0, max(max(iterations), max(iterations_qp))*1.1])
savefigpdf(f, "ex2_iter_comp_large", 2);




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





