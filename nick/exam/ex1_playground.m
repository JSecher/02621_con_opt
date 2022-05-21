H = [5.0000,1.8600,1.2400,1.4800 -0.4600;
    1.8600, 3.0000,  0.4400,  1.1200,  0.5200;
    1.2400,  0.4400,  3.8000,  1.5600, -0.5400;
    1.4800,  1.1200,  1.5600,  7.2000, -1.1200;
    -0.4600,  0.5200, -0.5400, -1.1200,  7.8000];

g = [-16.1000, -8.5000, -15.7000, -10.0200, -18.6800];

A = [16.1000, 1.0000;
      8.5000, 1.0000;
     15.7000, 1.0000;
     10.0200, 1.0000;
     18.6800, 1.0000];

b = [10, 1];

N = 100000;
total_t1 = 0;
total_t2 = 0;

for i=1:N
    [x1, lam1, t1]= EqualityQPSolverNullSpace(H,g,A,b);
    total_t1 = total_t1 + t1(2);
end

for i=1:N
    [x2, lam2, t2]= EqualityQPSolverNullSpace(H,g,A,b);
    total_t2 = total_t2 + t2(2);
end


disp("time1")
disp(total_t1)
disp("time2")
disp(total_t2)


%%
[n,m] = size(A);

% Factorize A matrix using QR decomposition
tstart = cputime; 
[Q,Rbar] = qr(A, 'vector');
time1 = cputime-tstart;

if ~ m == size(Rbar, 2)
    error("A does not have full column rank")
end

% Solve for x and lagrange multipliers    
Q1 = Q(:,1:m);      % Q Range
Q2 = Q(:,m+1:n);    % Q Null
R = Rbar(1:m,1:m);
x_Y = R' \ b;       % Solve R' x_Y = b

% Solve: (Q2' H Q2)x_Z = −Q′2(HQ1 x_Y + g) i.e. solve for x_Z
Q2T = Q2';

%% 

N = 10000000;

total_t1 = 0;
total_t2 = 0;

tstart = cputime; 
for i=1:N
    (-Q2T * H*Q1*x_Y - Q2T * g);
end
time1 = cputime-tstart;

tstart = cputime; 
for i=1:N
    (-Q2T * (H*Q1*x_Y+g));
end
time2 = cputime-tstart;
disp("time1")
disp(time1)
disp("time2")
disp(time2)
%%

disp("time1")
disp(total_t1)
disp("time2")
disp(total_t2)


%% 
% Two axises 

% Plot total time used for for each method.
f = figure;
a1 = axes;
hold off
for i=1:length(solvers)
    plot(param_vals, times(i,:))
    hold on
end
legend(solvers, 'Location','northwest')

a2 = copyobj(a1,f);  % Create a new axis based on the existing axis
set(a2,'Color','none')  % Make the new axis transparent
set(a2,'Ytick',[])  % Remove the YTicks (since they are redundant)
set(a2,'XAxisLocation','top') % Display the XTicks on the top side
% set(a2,'XTickLabel', con_vals)
% Label your plot appropriately
ylabel(a1,"time [s]")
xlabel(a1,"n, number of parameters")
xlabel(a2,"m, number of constraints")