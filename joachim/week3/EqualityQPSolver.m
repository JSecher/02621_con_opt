function [x, lambda] = EqualityQPSolver(H, g, A, b)
    %EQUALITYQPSOLVER Function for solution of equality constrained convex quadratic programs.
    %   Detailed explanation goes here
    
    % If no constraints, just solve for x variabels
    if isempty(A) && isempty(b)
            x = H\-g;
            lambda = [];
    end
    
    [n, m] = size(H);
    KKT_mat = [H, -A; -A', zeros(length(b))];
    KKT_rhs = [-g; -b];

    res = KKT_mat \ KKT_rhs;
    x = res(1:m);
    lambda = res(m+1:end);

end

