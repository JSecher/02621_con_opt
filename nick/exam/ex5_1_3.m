%% Clean up
clear
close all

%% Exercise 5.1.3
%Solve for return of 12 and minimum risk
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

x = quadprog(H, [], Aineq, bineq, Aeq, beq)

port_risk = x'*cov*x