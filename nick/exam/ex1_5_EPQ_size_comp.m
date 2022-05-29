%% Problem 1.5 Driver 

%% Clean up
clear
close all

%% Driver for Exam Problem 1.5
% Investigation of size dependent problems
seed = 42;
Ntests = 20;
min_param = 100;
max_param = 2000; % 1000

run_param_dependent = true;
run_con_dependent = true;
run_recycle = true;

%% Run basic tests for all solvers
if run_param_dependent
    solvers = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', "quadprog"};
    % xs = zeros(length(solvers), Ntests, length(g));
    % lams = zeros(length(solvers), Ntests, length(g));
    err = zeros(length(solvers), Ntests);
    times = zeros(length(solvers), Ntests);
    times_factor= zeros(length(solvers), Ntests);
    
    param_vals = round(linspace(min_param, max_param, Ntests), -1); % round to nearest 10
    
    rng(seed)
    fprintf("Running %d tests, with Paramters from %d to %d, constraints %d to %d\n", Ntests, min_param, max_param, min_param, max_param)
    for j=1:Ntests
        fprintf("Iter %03d : n = %d \t m = %d\n", j, param_vals(j),param_vals(j))
        [H,g,A,b,x_true,lam_true] = randomEQP(param_vals(j),param_vals(j));
    
        for i=1:length(solvers)
            solver = solvers{i};
            tstart = cputime;
            if solver == "quadprog"
                x = quadprog(H,g,[],[],A',b, [], [], [], optimset('Display', 'off')); % Calling quadprog for true solution
            else
                [x, ~, times_factor(i,j)] = EqualityQPSolver(H, g, A, b, solver);
            end
            times(i,j) = cputime - tstart;
        
            % Compute error compared to true sol
            err(i,j) = norm(x - x_true);
        end
    end
    
    
    %% Plot basic run 
    
    % Plot total time used for for each method.
    f = figure("Name","Time Comparison");
    hold off
    grid on
    for i=1:length(solvers)
        plot(param_vals, times(i,:))
        hold on
    end
    legend(solvers, 'Location','northwest')
    ylabel("time [s]")
    xlabel("Number of parameters and constraints")
    savefigpdf(f, "ex1_5_time_comparison",1);
    
    % Plot total time used for for each method.
    f = figure("Name","Time Comparison, sparse removed");
    hold off
    grid on
    select_id = [1,3,5,6,7];
    for i=select_id
        plot(param_vals, times(i,:))
        hold on
    end
    legend(solvers(select_id), 'Location','northwest')
    ylabel("time [s]")
    xlabel("Number of parameters and constraints")
    savefigpdf(f, "ex1_5_time_comparison_no_sparse",1);
    
    % Plot total time used for for each method.
    f = figure("Name","Time Comparison, sparse removed, semilogy");
    hold off
    grid on
    select_id = [1,3,5,6,7];
    for i=select_id
        semilogy(param_vals, times(i,:))
        hold on
    end
    legend(solvers(select_id), 'Location','northwest')
    ylabel("time [s]")
    xlabel("Number of parameters and constraints")
    savefigpdf(f, "ex1_5_time_comparison_no_sparse_semilogy",1);
    
    % Plot total time used for for each method.
    f = figure("Name","Time Comparison, semilogy");
    hold off
    grid on
    for i=1:length(solvers)
        semilogy(param_vals, times(i,:))
        hold on
    end
    legend(solvers, 'Location','northwest')
    ylabel("time [s]")
    xlabel("Number of parameters and constraints")
    savefigpdf(f, "ex1_5_time_comparison_semilogy",1);
    
    % Plot total time used for for each method.
    f = figure("Name","Time Factorization Comparison");
    hold off
    grid on
    for i=1:length(solvers)
        plot(param_vals, times_factor(i,:))
        hold on
    end
    legend(solvers, 'Location','northwest')
    ylabel("time [s]")
    xlabel("Number of parameters and constraints")
    savefigpdf(f, "ex1_5_time_factorization_comparison",1);
    
    % Plot total error for each method.
    f = figure("Name","Error Comparison");
    hold off
    grid on
    for i=1:length(solvers)
        plot(param_vals, err(i,:));
        hold on
    end
    legend(solvers, 'Location','northwest')
    ylabel("MSE")
    xlabel("Number of parameters and constraints")
    savefigpdf(f, "ex1_5_err_comparison",1);
end

%% Growing constraints tests
if run_con_dependent
    % Run tests for all solvers
    solvers = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace'};
    con_err = zeros(length(solvers), Ntests);
    con_err_lam = zeros(length(solvers), Ntests);
    con_times = zeros(length(solvers), Ntests);
    con_times_factor= zeros(length(solvers), Ntests);
    
    con_vals = round(linspace(min_param, max_param, Ntests)); % round to nearest 10
    
    rng(seed)
    fprintf("Running %d tests, with fixed number of variables; %d, number of constraints fom %d to %d\n", Ntests, max_param, min_param, max_param)
    for j=1:Ntests
        fprintf("Iter %03d : n = %d \t m = %d\n", j, max_param,con_vals(j))
        [H,g,A,b,x_true,lam_true] = randomEQP(max_param,con_vals(j));
    
        for i=1:length(solvers)
            solver = solvers{i};
            tstart = cputime;
            [x, lam, con_times_factor(i,j)] = EqualityQPSolver(H, g, A, b, solver);
            con_times(i,j) = cputime - tstart;
        
            % Compute error compared to true sol
            con_err(i,j) = norm(x - x_true);
            con_err_lam(i,j) = norm(lam - lam_true);
        end
    end
    
    
    % Plot total time used for for each method.
    f = figure("Name","Growing constraints all");
    hold off
    grid on
    for i=1:length(solvers)
        plot(con_vals, con_times(i,:))
        hold on
    end
    legend(solvers, 'Location','northwest')
    ylabel("time [s]")
    xlabel("Number of constraints")
    savefigpdf(f, "ex1_5_constraint_time_comparison",1);

    % Plot total time used for for each method.
    f = figure("Name","Growing constraints Factorization all");
    hold off
    grid on
    for i=1:length(solvers)
        plot(con_vals, con_times_factor(i,:))
        hold on
    end
    legend(solvers, 'Location','northwest')
    ylabel("time [s]")
    xlabel("Number of constraints")
    savefigpdf(f, "ex1_5_constraint_time_factor_comparison",1);
    
    
    % Only nullspace and rangespace
    f = figure("Name","Growing constraints Select");
    hold off
    grid on
    solver_select = [5, 6];
    for i=solver_select
        plot(param_vals, times(i,:))
        hold on
    end
    legend(solvers(solver_select), 'Location','northwest')
    ylabel("time [s]")
    xlabel("Number of constraints")
    savefigpdf(f, "ex1_5_constraint_time_comparison_select",1);

end

%%%%%%%%%%%%%%%
%% Recycle problem
%%%%%%%%%%%%%%%
%% Do it all again for recycle problem
if run_recycle
    uhat = 0.2;
    d0 = 1;
    %% Run basic tests for all solvers recycle
    
    solvers = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace', "quadprog"};
    times = zeros(length(solvers), Ntests);
    times_factor= zeros(length(solvers), Ntests);
    
    param_vals = round(linspace(min_param, max_param, Ntests), -1); % round to nearest 10
    
    rng(seed)
    fprintf("Recycle: Running %d tests, with Paramters from %d to %d\n", Ntests, min_param, max_param);
    fprintf("\t uhat= %.3f, \t d0 = %.3f\n", uhat, d0);
    for j=1:Ntests
        fprintf("Iter %03d : n = %d \t m = %d\n", j, param_vals(j),param_vals(j))
        [H,g,A,b] = recycleEQP(param_vals(j), uhat, d0);
    
        for i=1:length(solvers)
            solver = solvers{i};
            tstart = cputime;
            if solver == "quadprog"
                x = quadprog(H,g,[],[],A',b, [], [], [], optimset('Display', 'off')); % Calling quadprog for true solution
            else
                [x, ~, times_factor(i,j)] = EqualityQPSolver(H, g, A, b, solver);
            end
            times(i,j) = cputime - tstart;
        end
    end
    
    
    %% Plot basic run Recyle
    
    % Plot total time used for for each method.
    f = figure("Name","Time Comparison");
    hold off
    grid on
    for i=1:length(solvers)
        plot(param_vals, times(i,:))
        hold on
    end
    legend(solvers, 'Location','northwest')
    ylabel("time [s]")
    xlabel("n")
    savefigpdf(f, "ex1_5_recycle_time_comparison",1);
    
    % Plot total time used for for each method.
    f = figure("Name","Time Comparison, sparse removed");
    hold off
    grid on
    select_id = [1,3,5,6,7];
    for i=select_id
        plot(param_vals, times(i,:))
        hold on
    end
    legend(solvers(select_id), 'Location','northwest')
    ylabel("time [s]")
    xlabel("n")
    savefigpdf(f, "ex1_5_recycle_time_comparison_no_sparse",1);
    
    % Plot total time used for for each method.
    f = figure("Name","Time Comparison, sparse removed, semilogy");
    hold off
    grid on
    select_id = [1,3,5,6,7];
    for i=select_id
        semilogy(param_vals, times(i,:))
        hold on
    end
    legend(solvers(select_id), 'Location','northwest')
    ylabel("time [s]")
    xlabel("n")
    savefigpdf(f, "ex1_5_recycle_time_comparison_no_sparse_semilogy",1);
    
    % Plot total time used for for each method.
    f = figure("Name","Time Comparison, semilogy");
    hold off
    grid on
    for i=1:length(solvers)
        semilogy(param_vals, times(i,:))
        hold on
    end
    legend(solvers, 'Location','northwest')
    ylabel("time [s]")
    xlabel("n")
    savefigpdf(f, "ex1_5_recycle_time_comparison_semilogy",1);
    
    % Plot total time used for for each method.
    f = figure("Name","Time Factorization Comparison");
    hold off
    grid on
    select_id = 1:(length(solvers)-1);
    for i=select_id
        plot(param_vals, times_factor(i,:))
        hold on
    end
    legend(solvers(select_id), 'Location','northwest')
    ylabel("time [s]")
    xlabel("n")
    savefigpdf(f, "ex1_5_recycle_time_factorization_comparison",1);

    % Plot total time used for for each method.
    f = figure("Name","Time Factorization Comparison, semilogy");
    hold off
    grid on
    select_id = 1:(length(solvers)-1);
    for i=select_id
        semilogy(param_vals, times_factor(i,:))
        hold on
    end
    legend(solvers(select_id), 'Location','northwest')
    ylabel("time [s]")
    xlabel("n")
    savefigpdf(f, "ex1_5_recycle_time_factorization_comparison_semilogy",1);
    
end





