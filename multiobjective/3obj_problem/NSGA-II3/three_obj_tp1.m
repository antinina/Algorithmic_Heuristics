function f = three_obj_tp1(X, alpha, beta)
%TP1_OBJ  Vectorized objective for the 3-objective test problem.
%   Input:
%       X : [N x n] population matrix (rows = candidates)
%   Output:
%       F : [N x 3] objective values [f1 f2 f3]

    % Ensure 2D

    if isvector(X), X = X(:).'; end

    x1 = X(:,1);
    x2 = X(:,2);
    n  = size(X,2);

    % h(x1,x2)
    h = 2 - x1.^2 - x2.^2;

    % g(x) = sum_{i=3..n} (10 + x_i^2 - 10 cos(4*pi*x_i))
    if n >= 3
        Xi = X(:,3:end);
        g  = sum( 10 + Xi.^2 - 10.*cos(4*pi*Xi), 2 );
    else
        g  = zeros(size(X,1),1);
    end

    % S(x1,x2)  (from the figure)
    % S(x1,x2) = alpha/(0.2 + x1) + beta*x1^8 + alpha/(0.2 + x2) + beta*x2^8
    S = alpha./(0.2 + x1) + beta.*(x1.^8) + alpha./(0.2 + x2) + beta.*(x2.^8);

    % Objectives
    f1 = x1;
    f2 = x2;
    f3 = h + g.*S;

    f = [f1, f2, f3]';
end