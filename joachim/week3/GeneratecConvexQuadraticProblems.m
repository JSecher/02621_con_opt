function [H,g,A,b,xtarget] = GeneratecConvexQuadraticProblems(maxcoef, n, m)
%CREATE Summary of this function goes here
%   Detailed explanation goes here
% maxcoef, max coefficent for f
% n, number of variabels
% m, number of constraints 

p = 1;
while p ~= 0
    d = (maxcoef+1)*rand(n,1); % The diagonal values
    t = triu(bsxfun(@min,d,d.').*rand(n,1),1); % The upper trianglar random values
    H = diag(d)+t+t.'; % Put them together in a symmetric matrix
    
    % Just to make sure it is postive definite
    [~,p] = chol(H);
end

xtarget = rand(n,1);
A = rand(n, m);

b = A' * xtarget;
g = -H * xtarget;


end

