%% Problem 5 Driver
%% Clean up
clear
close all

%% Section 5.3
% Set up market
returns = [16.1, 8.5, 15.7, 10.02, 18.68]; 
R_val = 12; 
Sigma = [2.5 0.93 0.62 0.74 -0.23;
         0.93 1.5 0.22 0.56 0.26;   
         0.62  0.22 1.9  0.78  -0.27;
         0.74 0.56 0.78 3.6 -0.56;
         -0.23 0.26 -0.27 -0.56 3.9];  

Aeq = [returns; [1,1,1,1,1]];
beq = [R_val; 1];
Aineq = -eye(5);
bineq = zeros(5,1); 

x_min = quadprog(Sigma, [], Aineq, bineq, Aeq, beq)

risk_12 = x_min'*Sigma*x_min

%% Section 5.4
% Ramping the return R for a 1000 values
options = optimset('Display', 'off');
n = 200;

R = linspace(8.5, 18.68, n);
risk_R = zeros(n,1);
port_R = zeros(n,5);

for i = 1:length(R)
beq = [R(i); 1];
x = quadprog(Sigma, [], Aineq, bineq, Aeq, beq,[],[],[],options); 
risk_R(i) = x'*Sigma*x;
port_R(i,:) = x;
end

%% Plotting efficient frontier
% illustrating the ramping 
f = figure("Name","Efficientcy frontier");
hold on

opt_return = R(risk_R==min(risk_R));
opt_risk = min(risk_R);

%Plot efficient frontier for each possible return
plot(R,risk_R,'r', 'LineWidth',2)
plot(opt_return, opt_risk ,'ko','MarkerFaceColor', 'black');
%%title('Efficient Frontier')
xlabel('Return')
ylabel('Risk [Var[R]]')
legend('Efficient frontier','Minimum risk portfolio', 'Location','northwest')
xlim([8.5, 18.68])
hold off

savefigpdf(f, "ex5_4_efficient_frontier",5);

%% Plotting the Portfolio with minimum variance as a function of return
f = figure("Name","Portfolio Distribution");
hold on
plot(R,port_R(:,1),"LineWidth",2);
plot(R,port_R(:,2),"LineWidth",2);
plot(R,port_R(:,3),"LineWidth",2);
plot(R,port_R(:,4),"LineWidth",2);
plot(R,port_R(:,5),"LineWidth",2);
legend('Security 1','Security 2','Security 3','Security 4','Security 5','Location','north')
xlabel('Return')
ylabel('Part of portfolio')
%%title('Portfolio distribution')
xlim([8.5, 18.68])
hold off
savefigpdf(f, "ex5_4_portfolio_distribution",5);

%% Pareto Optimal Point
beq = [opt_return; 1];
x_min = quadprog(Sigma, [], Aineq, bineq, Aeq, beq)
risk_min = x_min'*Sigma*x_min

%% Section 5.8-5.11
% Set up of extended financial market.
returns_ext = [16.1, 8.5, 15.7, 10.02, 18.68, 0]; 
R_val = 14; 
Sigma_ext = [2.5 0.93 0.62 0.74 -0.23 0;
         0.93 1.5 0.22 0.56 0.26 0;
         0.62  0.22 1.9  0.78  -0.27 0;
         0.74 0.56 0.78 3.6 -0.56 0;
         -0.23 0.26 -0.27 -0.56 3.9 0;
         0 0 0 0 0 0];  

Aeq_ext = [returns_ext; [1,1,1,1,1,1]];
beq = [R_val; 1];
Aineq_ext = -eye(6);
bineq_ext = zeros(6,1); 

x_min_ext = quadprog(Sigma_ext, [], Aineq_ext, bineq_ext, Aeq_ext, beq)

risk_14_ext = x_min_ext'*Sigma_ext*x_min_ext


%% Ramping R in the extended financial market
n = 200;

R_ext = linspace(0, 18.68, n);
risk_R_ext = zeros(n,1);
port_R_ext = zeros(n,6);

for i = 1:length(R_ext)
beq = [R_ext(i); 1];
x_ext = quadprog(Sigma_ext, [], Aineq_ext, bineq_ext, Aeq_ext, beq); 
risk_R_ext(i) = x_ext'*Sigma_ext*x_ext;
port_R_ext(i,:) = x_ext;
end

%% Plotting efficient frontier for extended market
f = figure("Name","Efficientcy frontier");
hold on

opt_return_ext = R_ext(risk_R_ext==min(risk_R_ext));
opt_risk_ext = min(risk_R_ext);

%Plot efficient frontier for each possible return
plot(R_ext,risk_R_ext,'r', 'LineWidth',2)
plot(opt_return_ext, opt_risk_ext ,'ko','MarkerFaceColor', 'black');
plot(returns_ext, diag(Sigma_ext), 'x', 'Color', 'black','MarkerSize',10);
%title('Efficient Frontier')
xlabel('Return')
ylabel('Risk [Var[R]]')
legend('Efficient frontier','Minimum risk portfolio','Security coordinates','Location','northwest')
xlim([0, 18.68])
hold off

savefigpdf(f, "ex5_8_efficient_frontier",5);

%% Plotting the Portfolio with minimum variance as a function of return
f = figure("Name","Portfolio Distribution");
hold on
plot(R_ext,port_R_ext(:,1),"LineWidth",2);
plot(R_ext,port_R_ext(:,2),"LineWidth",2);
plot(R_ext,port_R_ext(:,3),"LineWidth",2);
plot(R_ext,port_R_ext(:,4),"LineWidth",2);
plot(R_ext,port_R_ext(:,5),"LineWidth",2);
plot(R_ext,port_R_ext(:,6),"LineWidth",2);
xline(14.0, "--r",'LineWidth',2)
legend('Security 1','Security 2','Security 3','Security 4','Security 5','Security f' ,'Location','north')
xlabel('Return')
ylabel('Part of portfolio')
%title('Portfolio distribution')
xlim([0, 18.68])
hold off
grid off
savefigpdf(f, "ex5_8_portfolio_distribution",5);

%% Based on Frontier, we find the minimum risk (Pareto point)
beq = [opt_return_ext; 1];
x_min_ext = quadprog(Sigma_ext, [], Aineq_ext, bineq_ext, Aeq_ext, beq)
risk_min_ext = x_min_ext'*Sigma_ext*x_min_ext

%% Plotting efficient frontiers
f = figure("Name","Efficientcy frontiers");
hold on

plot(R,risk_R,'b', "LineWidth", 2)
plot(R_ext,risk_R_ext,'r', "LineWidth", 2)
plot(R_val, 0.7176, 'o', 'Color','blue','MarkerSize',10, 'Markerfacecolor', 'b');
plot(R_val, risk_14_ext, 'o', 'Color','red','MarkerSize',10, 'MarkerFaceColor','r');
%title('Efficient Frontier')
xlabel('Return')
ylabel('Risk [Var[R]]')
legend('Old Market','New Market', 'Location','northwest')
xlim([8.5, 18.68])

hold off

savefigpdf(f, "ex5_11_efficient_frontiers",5);

%% Section 5.5-5.7,
% Below we set up a Bi-criterion optimization program for the market.
f = -returns;

Aeq = [1,1,1,1,1];
beq = 1;

Aineq = -eye(5);
bineq = zeros(5,1); 
alpha_trials = 200;
x = zeros(alpha_trials,5,3);

alphas = linspace(0.001, 1, alpha_trials)
port_risk= zeros(alpha_trials,1,3);
port_return= zeros(alpha_trials,1,3);

x0 = zeros(5,1);
s0 = ones(2*5,1);
y0 = ones(length(beq),1);
z0 = ones(2*5,1);

quadtime = zeros(alpha_trials,1);
qpsolvertime = zeros(alpha_trials,1);
Nsamples = 5;

for i = 1:alpha_trials
    % Solving QP problem with buidin solver and our solver:
    tic
    for sample=1:Nsamples
    x(i,:,1) = quadprog( alphas(i).*Sigma, (1-alphas(i)).*f', Aineq, bineq, Aeq, beq,[],[],[],options); 
    end
    quadtime(i) = toc/Nsamples;
    tic
    for sample=1:Nsamples 
    x(i,:,2) = QPSolverMehrotraInteriorPoint(alphas(i).*Sigma, (1-alphas(i)).*f',Aeq',beq,zeros(5,1),ones(5,1),x0,y0,z0,s0);
    end
    qpsolvertime(i) = toc/Nsamples;

    port_risk(i,2) = x(i,:,2)*Sigma*x(i,:,2)';
    port_return(i,2) = -f*x(i,:,2)';
end


% MSE in comparrison to quadproc.
f = figure("Name","MSE compared to quadprog");
hold on
grid on
err = sum(sqrt((x(:,:,1)-x(:,:,2)).^2),2)
plot(alphas, err,'LineWidth',2);
xlabel('\alpha')
ylabel('MSE compared to quadprog')
hold off
%savefigpdf(f, "ex5_5_MSE_alpha",5);

f = figure("Name","Log10 Mean Squared Error compared with quadprog");
err = sum(sqrt((x(:,:,1)-x(:,:,2)).^2),2)
semilogy(alphas, err,'LineWidth',2);
hold on
grid on
xlabel('\alpha')
ylabel('MSE compared to quadprog')
hold off
%savefigpdf(f, "ex5_5_Log_MSE_alpha",5);

% Time - efficiency
f = figure("Name","Time comparrison");
hold on
grid on
plot(alphas, qpsolvertime,'LineWidth',2);
plot(alphas, quadtime,'LineWidth',2);
legend("Own Solver", "QuadProg")
%%title(sprintf("Mean over %d samples", Nsamples))
xlabel('\alpha')
ylabel('time [s]')
hold off
%savefigpdf(f, "ex5_5_time_comparisson",5);


%Efficient Frontier
f = figure("Name","Bi-Criterion");
hold on
plot(port_return(:,2),port_risk(:,2),'r', 'LineWidth',2)
xlim([8.5, 18.68])
%%title('Bi-Criterion sampling \alpha')
xlabel('Return')
ylabel('Risk [Var[R]]')
hold off
%savefigpdf(f, "ex5_5_bi_criterion_for_alpha",5);