%% Clean up
clear
close all

%% Exercise 5.4,
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
%% Plotting the plotting the efficient frontier
f = figure('Name','b vals Relative Error solver');
hold on
opt_return1 = Rs(pf_risks==min(pf_risks));
opt_risk1 = min(pf_risks);

effRs = Rs(find(Rs == opt_return1):end);
effRisk = pf_risks(find(Rs == opt_return1):end);

%Plot efficient frontier for each possible return
plot(Rs, pf_risks, returns, diag(cov),'ko')
grid on
p1 = plot(opt_return1,opt_risk1 ,'ko','MarkerFaceColor', 'g');
p2 = plot(effRs,effRisk , 'r');
title('The Efficient Frontier')
legend([p2,p1],{'The efficient frontier','Minimum efficient return'}, 'Location', 'northwest')
xlabel('Return [%]')
ylabel('Risk [Var]')
hold off
savefigpdf(f, "ex5_1_the_efficient_frontier", 5); 
%% Plotting the optimal portfolio composition as a function of return
figure
hold on 
grid on
p3 = plot(Rs,optPf(:,1));
p4 = plot(Rs,optPf(:,2));
p5 = plot(Rs,optPf(:,3));
p6 = plot(Rs,optPf(:,4));
p7 = plot(Rs,optPf(:,5));
legend([p3,p4,p5,p6,p7], {'Security 1','Security 2','Security 3','Security 4','Security 5'})
xlabel('Return [%]')
ylabel('Percentage of the portfolio [%]')
title('Optimal Portfolio Composition')
savefigpdf(f, "exe5_1_portfolio_coposisitions", 5);

%% Find the portfolio with the smallest possible risk.
beq = [opt_return1; b2];
x_min = quadprog(H, [], Aineq, bineq, Aeq, beq,[],[],[],options)
        
pf_risk_min = x'*cov*x

