%% Clean up
clear
clc
close all
%% Calculate what is needed from old market

% Computing the efficient frontier
returns = [16.1, 8.5, 15.7, 10.02, 18.68]; 
cov = [2.5 .93 .62 .74 -.23;
              0.93 1.5 0.22 0.56 0.26;
              .62  .22 1.9  .78  -0.27;
              .74 .56 .78 3.6 -0.56;
              -0.23 0.26 -0.27 -0.56 3.9];  
R = 12; 
H = cov;

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

%% Setting up our financial markets
returns = [16.1, 8.5, 15.7, 10.02, 18.68]; 
cov = [2.5 .93 .62 .74 -.23;
              0.93 1.5 0.22 0.56 0.26;
              .62  .22 1.9  .78  -0.27;
              .74 .56 .78 3.6 -0.56;
              -0.23 0.26 -0.27 -0.56 3.9];  
          
returns = [returns, 0];

cov = [cov, zeros(5,1)];
cov = [cov; zeros(1,6)];

H = cov;
f = []; 

A1 = returns;
b1 = min(A1);

A2 = [1,1,1,1,1,1];
b2 = 1;

Aeq = [A1; A2];
beq = [b1; b2];

Aineq = -eye(6);
bineq = zeros(6,1); 

%% Exercise 5.3.2 - 5.3.4,

% Computing the efficient frontier with a risk-free asset 
%Rs_free = min(returns):0.01:max(returns);
Rs_free = linspace(min(returns),max(returns),1500);
port_risks_free = zeros(length(Rs_free),1);
optimalPorts = zeros(length(Rs),6);

for i = 1:length(Rs_free)
    beq = [Rs_free(i); b2];
    x = quadprog( H, f, Aineq, bineq, Aeq, beq,[],[],[],options); 
    port_risks_free(i) = x'*cov*x;
    optimalPorts(i,:) = x;
end

%Plotting whe Eff Frontier
fig1 =  figure('Name','The Efficient frontier, risk free asset');
grid on
hold on
%plot(Rs_free,port_risks_free,'b', returns, diag(cov),'ro')
plot(Rs_free,port_risks_free,'b')
title('The Efficient Frontier, risk free asset')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off
savefigpdf(fig1, "ex5_risk_free_efficient_frontier", 5);


lower_risk_pareto = pf_risks(min(find(abs(round(port_risks_free-pf_risks,2))==0)));
lower_return_pareto = Rs(min(find(abs(round(port_risks_free-pf_risks,2))==0)));


fig2 =  figure('Name','The Efficient frontier, risk free asset');
grid on
hold on
h1 = plot(Rs_free,port_risks_free,'b');
%scatter( returns, diag(cov),'ro')
h3 = plot(Rs,pf_risks,'color',[0 0.5 0]);
plot(opt_return1, opt_risk1,'ko','MarkerFaceColor', 'k', 'MarkerSize',3);
plot(lower_return_pareto,lower_risk_pareto,'ko','MarkerFaceColor', 'k', 'MarkerSize',3);
h2 = plot(effRs,effRisk , 'r');
hold off
title('The Efficient Frontiers, comparison')
xlabel('Return [%]')
ylabel('Risk [Var]')
legend([h1,h2, h3],{'With risk free asset','Without risk free asset and pareto optimal','Without risk free asset and not pareto optimal'}, 'Location', 'northwest')
savefigpdf(fig2, "ex5_frontier_market_comparison", 5);

%% Exercise, Portfolio compositions
% Plotting the optimal portfolio composition in new market
fig3 =  figure('Name','PF composition');
h1 = plot(Rs_free,optimalPorts(:,1));
hold on
h2 = plot(Rs_free,optimalPorts(:,2));
h3 = plot(Rs_free,optimalPorts(:,3));
h4 = plot(Rs_free,optimalPorts(:,4));
h5 = plot(Rs_free,optimalPorts(:,5));
h6 = plot(Rs_free,optimalPorts(:,6));
legend([h1,h2,h3,h4,h5, h6], {'Security 1','Security 2','Security 3','Security 4','Security 5', 'Risk free asset'})
xlabel('Return [%]')
ylabel('Percentage of the portfolio [%]')
title('Portfolio composition')
savefigpdf(fig3, "ex5_new_market_composition", 5);


%% Finding the optimum for R=14

beq = [14; b2];
x = quadprog( H, f, Aineq, bineq, Aeq, beq); 
opt_risk2 = x'*cov*x;

%% Plot the location R=14,  with and without a risk-free asset,
fig4 =  figure('Name','R=14 location with and without');
grid on
hold on
opt_return2 = 14;
%plot(Rs_free,port_risks_free,'b', returns, diag(cov),'mo')
plot(Rs_free,port_risks_free,'b')
plot(Rs,pf_risks,'color',[0 0.5 0])
h2 = plot(effRs,effRisk , 'r');
p1 = plot(opt_return2, opt_risk2,'ko','MarkerFaceColor', 'g', 'MarkerSize',4);
%plot(opt_return1, opt_risk1,'k|','MarkerFaceColor', 'k', 'MarkerSize',8);
plot(opt_return1, opt_risk1, 'k*');
[~,idx] = (min(abs(Rs - 14)));
p2 = plot(Rs(idx), pf_risks(idx),'ko','MarkerFaceColor', 'k', 'MarkerSize',4);
title('Comparison of Efficient Frontiers')
xlabel('Return [%]')
ylabel('Risk [Var]')
xlim([10 19])
ylim([0,1.5])
legend([p1,p2],{'E[R]=14 with risk free asset','E[R]=14 without risk free asset'},'Location','northwest') 
hold off
savefigpdf(fig4, "ex5_R14_with_without", 5);