%% Exercise 5.2.1-5.2.3,

%% Clean up
clear
close all

%% Setting up our financial market
returns = [16.1, 8.5, 15.7, 10.02, 18.68]; 
cov = [2.5 .93 .62 .74 -.23;
              0.93 1.5 0.22 0.56 0.26;
              .62  .22 1.9  .78  -0.27;
              .74 .56 .78 3.6 -0.56;
              -0.23 0.26 -0.27 -0.56 3.9];  
R = 12; 
H = cov;
f = -returns;
A1 = returns;
b1 = R;

A2 = [1,1,1,1,1];
b2 = 1;

Aeq = [A1; A2];
beq = [b1; b2];

Aineq = -eye(5);
bineq = zeros(5,1); 

Rs = linspace(min(returns),max(returns),1500);
pf_risks = zeros(length(Rs),1);
options = optimset('Display', 'off');
optPf = zeros(length(Rs),5);

%Compute the efficient frontier for each possible return
for i = 1:length(Rs)
beq = [Rs(i); b2];
x = quadprog( H, [], Aineq, bineq, Aeq, beq,[],[],[],options); 
pf_risks(i) = x'*cov*x;
optPf(i,:) = x;
end

opt_return1 = Rs(pf_risks==min(pf_risks));
opt_risk1 = min(pf_risks);

effRs = Rs(find(Rs == opt_return1):end);
effRisk = pf_risks(find(Rs == opt_return1):end);

%% Setup for bi-criterion
H = cov;
f = -returns;

Aeq = [1,1,1,1,1];
beq = 1;

Aineq = -eye(5);
bineq = zeros(5,1);  

%% Plot the efficient frontiers portfolios when 

%allowing or disallowing shorting, solving the Bicriterion
n_runs = 1000;
x_nonshort = zeros(n_runs,5,3);
x_short = zeros(n_runs,5,3);


alphas = linspace(0.001,1,n_runs);
port_risk_nonshort = zeros(n_runs,1,3);
port_risk_short = zeros(n_runs,1,3);


port_return_nonshort = zeros(n_runs,1,3);
port_return_short = zeros(n_runs,1,3);

x0 = zeros(5,1);
s0 = ones(2*5,1);
y0 = ones(length(beq),1);
z0 = ones(2*5,1);
cvx_solve = true;

n = 5;
for i = 1:n_runs 
    if mod(i,n_runs/10)==0 %Keep track how far we are :)
        disp(i)
    end
    
    % Without shorting, i.e. an IQP problem
    x_nonshort(i,:,1) = quadprog(alphas(i).*H, (1-alphas(i)).*f', Aineq, bineq, Aeq, beq,[],[],[]); 
    x_nonshort(i,:,2) = QPSolverMehrotraInteriorPoint(alphas(i).*H, (1-alphas(i)).*f',Aeq',beq,zeros(5,1),ones(5,1),x0,y0,z0,s0);
    if cvx_solve
        cvx_begin quiet
            %cvx_precision low
            variable x(n)
            minimize( 1/2 * x' *  (alphas(i).*H) *x + ((1-alphas(i)).*f)*x)
            subject to
                Aeq * x == beq
                zeros(5,1) <= x
                x <= ones(5,1)
        cvx_end
        x_nonshort(i,:,3) = x; 
    end
    
    % With shorting, i.e. an EQP problem
    x_short(i,:,1) = quadprog(alphas(i).*H, (1-alphas(i)).*f', [], [], Aeq, beq,[],[],[]); 
    x_short(i,:,2) = EqualityQPSolver(alphas(i).*H, (1-alphas(i)).*f', Aeq', beq, "rangespace");
                   
    if cvx_solve
        cvx_begin quiet
            %cvx_precision low
            variable x(n)
            minimize( 1/2 * x' *  (alphas(i).*H) *x + ((1-alphas(i)).*f)*x)
            subject to
                Aeq * x == beq
        cvx_end
        x_short(i,:,3) = x; 
    end
    
    port_risk_nonshort(i,2) = x_nonshort(i,:,2)*cov*x_nonshort(i,:,2)';
    port_risk_short(i,2) = x_short(i,:,2)*cov*x_short(i,:,2)';
    
    port_return_nonshort(i,2) = -f*x_nonshort(i,:,2)';
    port_return_short(i,2) = -f*x_short(i,:,2)';
end

%Plot the efficient frontier of the bi-criterion problem when shorting is
%allowed
fig1 = figure('Name','Bi-Criterion, shorting allowed');
hold on
plot(port_return_short(:,2),port_risk_short(:,2),'b', returns, diag(cov),'ro')
title('The Bi-Criterion when shorting is allowed')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off
savefigpdf(fig1, "ex5_risk-return-shorting", 5);

%And the efficient frontier of the Bi-Criterion problem when shorting is
%not allowed
fig2 = figure('Name','Bi-Criterion, no shorting allowed');
hold on
plot(port_return_nonshort(:,2),port_risk_nonshort(:,2),'b', returns, diag(cov),'ro')
title('The Bi-Criterion when shorting is not allowed')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off
savefigpdf(fig2, "ex5_risk-return-no-shorting", 5);

% Plot comparing the solutions to the Bi-Criterion with no shorting and
% optimizing for risk given a return.
fig3= figure('Name','Risk-Return Optim');
hold on
h1 = plot(Rs,pf_risks,'r');
plot(returns, diag(cov),'ro')
h2 = plot(port_return_nonshort(:,2),port_risk_nonshort(:,2),'b');
plot(opt_return1,opt_risk1 ,'k*','MarkerSize', 14)
title('Combined static return and bi-criterion')
xlabel('Return [%]')
ylabel('Risk [Var]')
legend([h1,h2],{'Static return','Bi-criterion'}, 'Location', 'northwest')
hold off
savefigpdf(fig3, "ex5_risk_return_static_return", 5);

%% Exercise 5.5-5.7, Solving for different alpha with and without short-selling
n_runs = 500;
x_short = zeros(n_runs,5,3);
x_nonshort = zeros(n_runs,5,3);


alphas = linspace(0.001,1,n_runs);
port_risk_short = zeros(n_runs,1,3);
port_risk_nonshort = zeros(n_runs,1,3);
port_return_short = zeros(n_runs,1,3);
port_return_nonshort = zeros(n_runs,1,3);

x0 = zeros(5,1);
s0 = ones(2*5,1);
y0 = ones(length(beq),1);
z0 = ones(2*5,1);
n = 5;

clear log
for i = 1:n_runs
    if mod(i,n_runs/10)==0
        disp(i)
    end
    % Without shorting, i.e. an IQP problem
    x_nonshort(i,:,1) = quadprog( alphas(i).*H, (1-alphas(i)).*f', Aineq, bineq, Aeq, beq,[],[],[],options); 
    x_nonshort(i,:,2) = QPSolverMehrotraInteriorPoint(alphas(i).*H, (1-alphas(i)).*f',Aeq',beq,zeros(5,1),ones(5,1),x0,y0,z0,s0);
    if cvx_solve
        cvx_begin quiet
            %cvx_precision low
            variable x(n)
            minimize( 1/2 * x' *  (alphas(i).*H) *x + ((1-alphas(i)).*f)*x)
            subject to
                Aeq * x == beq
                zeros(5,1) <= x
                x <= ones(5,1)
        cvx_end
        x_nonshort(i,:,3) = x; 
    end
    
    % With shorting, i.e. an EQP problem
    x_short(i,:,1) = quadprog( alphas(i).*H, (1-alphas(i)).*f', [], [], Aeq, beq,[],[],[],options); 
    x_short(i,:,2) = EqualityQPSolver(alphas(i).*H, (1-alphas(i)).*f',Aeq', beq, "rangespace");
    if cvx_solve
        cvx_begin quiet
            %cvx_precision low
            variable x(n)
            minimize( 1/2 * x' *  (alphas(i).*H) *x + ((1-alphas(i)).*f)*x)
            subject to
                Aeq * x == beq
        cvx_end
        x_short(i,:,3) = x; 
    end
    
    port_risk_short(i,2) = x_short(i,:,2)*cov*x_short(i,:,2)';
    port_risk_nonshort(i,2) = x_nonshort(i,:,2)*cov*x_nonshort(i,:,2)';
    
    port_return_short(i,2) = -f*x_short(i,:,2)';
    port_return_nonshort(i,2) = -f*x_nonshort(i,:,2)';
end


%Plot the MSE for our method, quadprog, and cvx, solving the EQP problem
fig1 = figure('Name','Comparison of solvers no shorting');
hold on
h1 = plot((1:n_runs)./n_runs,log10(sum(sqrt((x_short(:,:,1)-x_short(:,:,2)).^2),2)));
hold on 
h2 = plot((1:n_runs)./n_runs,log10(sum(sqrt((x_short(:,:,3)-x_short(:,:,2)).^2),2)));
h3 = plot((1:n_runs)./n_runs,log10(sum(sqrt((x_short(:,:,3)-x_short(:,:,1)).^2),2)));
title('Comparison of found portfolios with shorting')
xlabel('alpha')
ylabel('Mean Squared Error, log10')
legend([h1,h2,h3],{'MSE log10, quadprog and our','MSE log10, CVX and our','MSE log10, CVX and quadprog'}, 'Location', 'northeast')
savefigpdf(fig1, "ex5_compare_shorting", 5);

%And when solving the full QP problem
fig2 = figure('Name','Comparison of solvers with shorting');
h1 = plot((1:n_runs)./n_runs,log10(sum(sqrt((x_nonshort(:,:,1)-x_nonshort(:,:,2)).^2),2)));
hold on 
h2 = plot((1:n_runs)./n_runs,log10(sum(sqrt((x_nonshort(:,:,3)-x_nonshort(:,:,2)).^2),2)));
h3 = plot((1:n_runs)./n_runs,log10(sum(sqrt((x_nonshort(:,:,3)-x_nonshort(:,:,1)).^2),2)));
title('Comparison of found portfolios with no shorting')
xlabel('alpha')
ylabel('Mean Squared Error, log10')
legend([h1,h2,h3],{'MSE log10, quadprog and our','MSE log10, CVX and our','MSE log10, CVX and quadprog'}, 'Location', 'northeast')
savefigpdf(fig2, "ex5_compare_no_shorting", 5);

