%% program using built in ga solver
nvars = 3; 

% Lower and upper bounds
lb = [-1*ones(1,nvars)];
ub = [ 1*ones(1,nvars)];

lb(1) = 0;   % first decision variable between 0 and 1
ub(1) = 1;

opts = optimoptions('gamultiobj', ...
    'PopulationSize', 50, ...
    'MaxGenerations', 100, ...
    'FunctionTolerance', 1e-8, ...
    'ParetoFraction', 0.7, ...
    'UseVectorized', false, ...
    'PlotFcn', {@gaplotpareto}, ...
    'Display', 'iter');


%% Run multiobjective GA
[x, fval] = gamultiobj(@TP1, nvars, [], [], [], [], lb, ub,[],opts);

%% plotting pareto front
figure;
plot(fval(:,1), fval(:,2), 'ro')
xlabel('Objective 1 (f1)')
ylabel('Objective 2 (f2)')
title('Pareto Front for TP1')
grid on
hold on

%% original efficient front

ff1 = linspace(-1,1,200);   
ff2 = 1 - ff1.^2;            
plot(ff1, ff2, 'k-', 'LineWidth',1, 'DisplayName','original efficient front');

hold off;
