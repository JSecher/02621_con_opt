%% Problem 1.4 Driver
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

% Generate values for b(1) in the range [8.5 18.68]
Ntests = 30;
bs = linspace(8.5, 18.68, Ntests);

%% Run tests for all solvers
solvers = {'QuadProg', 'LUdense', 'LUsparse', 'LDLdense', 'LDLsparse', 'RangeSpace', 'NullSpace'};
xs = zeros(length(solvers), Ntests, length(g));
err = zeros(length(solvers), Ntests);

for i=1:length(solvers)
    solver = solvers{i};
    fprintf("Using solver: %s\n", solver)
    for j=1:Ntests
        b(1) = bs(j);
        if i == 1   % First solver should be quadprog, this is used for error
            xs(i,j,:) = quadprog(H,g,[],[],A',b, [], [], [], optimset('Display', 'off')); % Calling quadprog for true solution
        else
            xs(i,j,:) = EqualityQPSolver(H, g, A, b, solver);
        end     
        % Compute error compared to quadprog
        err(i,j) = mean(sqrt((squeeze(xs(i,j,:) - xs(1,j,:)).^2)));
    end
end


%% Plot Error compared to QuadProg
f = figure('Name','Relative Error of solvers');
f.Position(3:4) = [800, 400];
for i=2:length(solvers)
    plot(bs, (err(i,:)), 'LineWidth',2)
    grid on
    hold on
end
legend(solvers(2:end))
xlabel("b(1)")
ylabel("MSE compared to quadprog")
savefigpdf(f, "ex1_4_error_comp_quadprog", 1);

%% Make Table

no_ids = 5;
solver_select = 2;
select_ids = round(linspace(1, Ntests, no_ids));
objective_true = zeros(no_ids,1);
objective = zeros(no_ids,1);

for i=1:length(select_ids)
    id = select_ids(i);
    x_quad = squeeze(xs(1,id,:));
    x_est = squeeze(xs(solver_select,id,:));
    objective_true(i) = 0.5*x_quad'*H*x_quad + g'*x_quad;
    objective(i) = 0.5*x_est'*H*x_est + g'*x_est;
end


data = objective_true';
data(2,:) = objective';
data(3,:) = (objective - objective_true)';
data(4,:) = nan(size(objective_true'));
data(5, :) = squeeze(err(solver_select,select_ids));
input.data = data;

% Set column labels (use empty string for no label):
input.tableColLabels = sprintfc('$\\beta$=%.4f',bs(select_ids));
% Set row labels (use empty string for no label):
input.tableRowLabels = {'$\phi_{\text{quadprog}}=$', ...
                        sprintf('$\\phi_{\\text{%s}}=$', solvers{solver_select}), ...
                        '$\Delta\phi=$', ...
                        '', ...
                        'MSE'};
% Set the row format of the data values 
input.dataFormat = {'%.4e'};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';
% Switch table borders on/off:
input.booktabs = 1;
% LaTex table caption:
input.tableCaption = sprintf('Comparison in ebjective value between \\textit{quadprog} and \\textit{%s}.', solvers{solver_select});
% LaTex table label:
input.tableLabel = 'ex1_4_objective_comp';
input.makeCompleteLatexDocument = 0;
input.dataNanString = '';
input.tablePlacement = 'ht';
% Now call the function to generate LaTex code:
latex = latexTable(input);
savelatexTable(latex, input.tableLabel, 1);





