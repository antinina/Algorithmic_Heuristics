function [GD, Delta] = ComputeMetrics(P, PFtrue)
% ComputeGDSpreed: Computes Generational Distance and Spread metrics
%
% Inputs:
%   P       - m x N matrix of obtained non-dominated solutions (each column = point)
%   PFtrue  - m x M matrix of true Pareto front points
%
% Outputs:
%   GD      - Generational Distance (convergence)
%   Delta   - Spread / diversity metric






%% 1. Generational Distance (GD)
    N = size(P,2);
    GDsum = 0;
    for i = 1:N
        pi = P(:,i);
        distances = sqrt(sum((PFtrue - pi).^2, 1)); % Euclidean distance to all true PF points
        d_i = min(distances);
        GDsum = GDsum + d_i^2;
    end
    GD = sqrt(GDsum / N);

    %% 2. Spread / Diversity (Δ)
    % Sort by first objective
    [~, idx] = sort(P(1,:));
    P = P(:, idx);

    % Distances between consecutive solutions
    d = sqrt(sum(diff(P,1,2).^2,1));
    dbar = mean(d);

    % Extreme points of true Pareto front
    [~, idxTrue] = sort(PFtrue(1,:));
    f_ext1 = PFtrue(:, idxTrue(1));
    f_ext2 = PFtrue(:, idxTrue(end));

    % Distances from obtained extremes to true extremes
    df = norm(P(:,1) - f_ext1);
    dl = norm(P(:,end) - f_ext2);

    % Spread formula
    Delta = (df + dl + sum(abs(d - dbar))) / (df + dl + (numel(d))*dbar + eps);

end
