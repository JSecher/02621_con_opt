%% Problem 3.4 - 3.6 Driver
%% Clean up
clear
close all

%% Test on toy problem for correctness

Nsamples = 1; % Run test multiple times for small problems
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
default_options = optimset('Display', 'off');

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

err = mean(sqrt((x - x_linprog).^2));
fprintf("Linprog Algo | Obj: %.5f \t Time: %.5f \t Iter: %3d\n", obj_linprog, time_linprog, iter_linprog)
fprintf("Own Algo     | Obj: %.5f \t Time: %.5f \t Iter: %3d \t MSE %.5f\n", obj, time_own, iter, err)
fprintf("Simplex      | Obj: %.5f \t Time: %.5f \t Iter: %3d\n", obj_simplex, time_simplex, iter_simplex)

fprintf("\n Thats nice...!, lets make some latex out of that, \n \n * *abrakadabra*\n \n *pafff*\n")

%% Make Table


data = [obj_linprog, obj];
data(2:6,:) = [x_linprog, x];
data(7,:) = [nan, err];
data(8,:) = [nan, nan];
data(9,:) = [iter_linprog, iter];
data(10,:) = [time_linprog, time_own];
input.data = data;

% Set column labels (use empty string for no label):
input.tableColLabels = {'Linprog IP', 'Own Solver'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'$\phi =$', ...
                        'x =', '', '', '', '', ...
                        'MSE =', ...
                        '', ...
                        'Iterations =','Time ='};
% Set the row format of the data values 
input.dataFormat = {'%.4f'};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'r';
% Switch table borders on/off:
input.booktabs = 1;
% LaTex table caption:
input.tableCaption = sprintf('Comparison between LinProg interior-point and Own LP interior-point alogorithm.');
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
data(2:6,:) = [x_linprog, x, x_simplex];
data(7,:) = [nan, err, nan];
data(8,:) = [nan, nan, nan];
data(9,:) = [iter_linprog, iter, iter_simplex];
data(10,:) = [time_linprog, time_own, time_simplex];
input.data = data;


% Set column labels (use empty string for no label):
input.tableColLabels = {'Linprog IP', 'Own Solver', 'Linprog Simplex'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'$\phi =$', ...
                        'x =', '', '', '', '', ...
                        'MSE =', ...
                        '', ...
                        'Iterations =','Time ='};
% Set the row format of the data values 
input.dataFormat = {'%.4f'};
% Column alignment ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'r';
% Switch table borders on/off:
input.booktabs = 1;
% LaTex table caption:
input.tableCaption = sprintf('Comparison between LinProg interior-point, Own LP interior-point and LinProg Simplex alogorithm.');
% LaTex table label:
input.tableLabel = 'ex3_objective_comp_simplex';
input.makeCompleteLatexDocument = 0;
input.dataNanString = '';
input.tablePlacement = 'ht';
% Now call the function to generate LaTex code:
latex = latexTable(input);
savelatexTable(latex, input.tableLabel, 3);
