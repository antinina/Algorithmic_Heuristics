% Three-objective test problem
% Minimizes: f1 = x1, f2 = x2, f3 = h(x1,x2) + g(x)*S(x1,x2)
% Bounds: 0<=x1,x2<=1 ; -1<=xi<=1 for i=3..n

clear;
clc;


n      = 3;         % 
alpha  = 0.75;
beta   = 10;

% Bounds 0<=x1,x2<=1 ; -1<=xi<=1 for i=3..n
% Lower and upper bounds
lb = [-1*ones(1,n)];
ub = [ 1*ones(1,n)];

lb(1) = 0;   % first variable between 0 and 1
ub(1) = 1;
lb(2) = 0;   % 2nd variable between 0 and 1
ub(2) = 1;


% lb = [-inf(1,n)];
% ub = [inf(1,n)];   % temp, then fix per-index
% lb(1:2) = [0 0];
% ub(1:2) = [1 1];
% lb(3:n) = -1;
% ub(3:n) =  1;

% ---------- Solver options ----------
opts = optimoptions('gamultiobj', ...
    'PopulationSize', 400, ...
    'MaxGenerations', 300, ...
    'FunctionTolerance', 1e-8, ...
    'ParetoFraction', 0.7, ...
    'UseVectorized', true, ...
    'PlotFcn', {@gaplotpareto}, ...
    'Display', 'iter');

% ---------- Solve ----------
fitness = @(X) three_obj_tp1(X, alpha, beta);  % vectorized objective
[xPF, fPF] = gamultiobj(fitness, n, [], [], [], [], lb, ub, [], opts);
%Mozda bolji prikaz
figure;
plot3(fPF(:,1), fPF(:,2), fPF(:,3), '.', 'MarkerSize', 14);
grid on; box on;
xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
title('Pareto Front approximation');
%% ovo za pareto koji je vec poznat
[FF1, FF2] = meshgrid(linspace(0,1,51), linspace(0,1,51));
FF3 = 2 - FF1.^2 - FF2.^2;

hold on;
mesh(FF1, FF2, FF3, ...
    'EdgeAlpha', 0.25, 'FaceAlpha', 0.15, 'FaceColor','interp');
legend('GA Pareto solutions','Analytical front');
%% Dots → Approximate Pareto-optimal solutions found by GA, Surface → True Pareto-optimal front.



if n > 2
    fprintf('Mean |x3..n| at GA solutions: %.3g\n', mean(abs(xPF(:,3:end)), 'all'));
end





