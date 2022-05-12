%% Problem 1.5 Driver 

%% Clean up
clear
close all

%% Driver for Exam Problem 1.5
% Investigation of size dependent problems


% Generate values for b(1) in the range [8.5 18.68]
Ntests = 30;
min_params = 10;
max_params = 1000;
min_con = min_params;
max_con = max_params;
seed = 42;


%% Run tests for all solvers
solvers = {'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace'};
% xs = zeros(length(solvers), Ntests, length(g));
% lams = zeros(length(solvers), Ntests, length(g));
err = zeros(length(solvers), Ntests);
err_lam = zeros(length(solvers), Ntests);
times = zeros(length(solvers), Ntests);
times_factor= zeros(length(solvers), Ntests);

param_vals = round(linspace(min_params, max_params, Ntests), -1); % round to nearest 10
con_vals = round(linspace(min_params, max_params, Ntests), -1); % round to nearest 10

rng(seed)
fprintf("Running %d tests, with Paramters from %d to %d, constraints %d to %d\n", Ntests, min_params, max_params, min_con, max_con)
for j=1:Ntests
    fprintf("Iter %03d : m = %d \t m = %d\n", j, param_vals(j),con_vals(j))
    [H,g,A,b,x_true,lam_true] = randomEQP(param_vals(j),con_vals(j));

    for i=1:length(solvers)
        solver = solvers{i};
        tstart = cputime;
        [x, lam, times_factor(i,j)] = EqualityQPSolver(H, g, A, b, solver);
        times(i,j) = cputime - tstart;
    
        % Compute error compared to true sol
        err(i,j) = norm(x - x_true);
        err_lam(i,j) = norm(lam - lam_true);
    end
end

%% Recycling problem
names = ["LDLdense", "LDLsparse", "LUdense", "LUsparse", "rangespace", "nullspace", "quadprog"];
tests = 20;
times = ones(7,tests)*tests;
ns = zeros(tests,1);
%We test every solver using the Recycling problem of different sizes n.
for i = 1:tests
    disp(i)
    n = i*(2000/tests);
    ns(i) = n;

    [H, g, A, b] = ProblemEQPRecycling(n, 0.2, 1);
    for k =1:1
        j=1;

        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LDLdense");
        times(j, i) = cputime-start;
        j = j+1;

        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LDLsparse");
        times(j, i) = cputime-start;
        j = j+1;
        %
        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LUdense");
        times(j, i) = cputime-start;
        j = j+1;
        %
        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "LUsparse");
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        %
        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "rangespace");
        times(j, i) = cputime-start;
        j = j+1;

        start = cputime;
        [x, lambda] = EqualityQPSolver(H,g,A,b, "nullspace");
        times(j, i) = cputime-start;
        j = j+1;

        options = optimset('Display', 'off');
        start = cputime;
        [x, lambda] = quadprog(H,g,[],[],A',b,[],[],[],options);
        times(j, i) = cputime-start;
        j = j+1;
        %
    end

    disp("n="+n)
end

%We plot the time for each method.
figure
hold off
for i=1:size(times,1)
    plot(ns, times(i,:))
    hold on
end
legend(names, 'Location','northwest')
xlabel("n")
ylabel("time [s]")
title("Growing size problem")


%% Factorization benchmark
tests = 10;
times = ones(5,tests)*100;
ns = zeros(tests,1);
%We bench LU vs LDL versus Cholesky of varying sizes
for i = 1:tests
    n = i*(5000/tests);
    ns(i) = n;

    [H,g,A,b,x,lambda] = generateRandomEQP(n,n);
    KKT = [H -A;-A', zeros(size(A,2), size(A,2))];
    for k =1:1
        j=1;
        
        start = cputime;
        x = lu(KKT,'vector');
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        x = lu(H,'vector');
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        x = ldl(KKT,'vector');
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        x = ldl(H,'vector');
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        x = chol(H);
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;


    end

    disp("n="+n)
end


%We then plot the time for each method and target.
hold off
for i=1:5%size(times,1)
    semilogy(ns, times(i,:))
    hold on
end
legend(["LU, KKT", "LU, H", "LDL, KKT", "LDL, H", "Cholesky, H"])
xlabel("n")
ylabel("time [s]")
title("Benchmark of factorizations")

%% m dependent
names = ["rangespace", "nullspace"];
tests = 10;
times = ones(2,tests)*100;
ms = zeros(tests,1);
top = 3000;
%To compare range space and null space, we solve random EQPs with varying
%number of constraints.
for i = 1:tests
    m = i*(top/tests);
    ms(i) = m;

    [H, g, A, b] = generateRandomEQP(top, m);
    for k =1:1
        j=1;


        start = cputime;
        [x, lambda, facTime_R] = EqualityQPSolver(H,g,A,b, "rangespace");
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

        start = cputime;
        [x, lambda, facTime_N] = EqualityQPSolver(H,g,A,b, "nullspace");
        times(j, i) = min(cputime-start, times(j,i));
        j = j+1;

    end

    disp("m="+m)
end

figure
hold off
for i=1:2
    plot(ms, times(i,:))
    hold on
end
plot(ms, times(1,:)-facTime_R)
plot(ms, times(2,:)-facTime_N)
xline(1950)
legend([names,'range space -Cholesky','null space -QR','Theoretical tipping point'], 'Location','northwest')
xlabel("m")
ylabel("time [s]")
title("Growing constraints problem, n=3000")

