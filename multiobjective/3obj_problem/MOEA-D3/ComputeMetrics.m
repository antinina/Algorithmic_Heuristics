function [GD, Delta] = ComputeMetrics(P, PFtrue)
% ComputeGDSpread3D: Computes GD and Spread metrics for any number of objectives
%
% Inputs:
%   P       - m x N matrix of obtained non-dominated solutions (columns = points)
%   PFtrue  - m x M matrix of true Pareto front points (columns = points)
%
% Outputs:
%   GD      - Generational Distance (lower is better)
%   Delta   - Spread / diversity metric (lower is better)

    %% 1. Generational Distance (GD)
    N = size(P,2);
    GDsum = 0;
    for i = 1:N
        pi = P(:,i);
        distances = sqrt(sum((PFtrue - pi).^2, 1));  % Euclidean distance to all true PF points
        d_i = min(distances);                         % nearest distance
        GDsum = GDsum + d_i^2;
    end
    GD = sqrt(GDsum / N);

    %% 2. Spread / Diversity (Delta)
    % Compute pairwise Euclidean distances between all solutions
    D = pdist2(P', P');   % N x N
    D_sorted = sort(D, 2); % sort each row
    d_nn = D_sorted(:,2);  % nearest neighbor distance for each point
    dbar = mean(d_nn);

    % Delta formula for m >= 3
    Delta = sum(abs(d_nn - dbar)) / (sum(d_nn) + eps);

end
