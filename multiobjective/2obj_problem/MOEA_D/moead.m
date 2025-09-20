clc;
clear;
close all;

%% Problem Definition

CostFunction = @(x) TP1(x);   % Cost Function

nVar = 3;                     % Number of Decision Variables (change as needed)
VarSize = [nVar 1];           % Decision Variables Matrix Size

% Variable Bounds
VarMin = [0; -1*ones(nVar-1,1)];   % Lower bounds: x1>=0, others >= -1
%VarMin = [0; 0*ones(nVar-1,1)];   % Lower bounds: x1>=0, others >= -1
VarMax = [1;  1*ones(nVar-1,1)];   % Upper bounds: x1<=1, others <= 1

nObj = numel(CostFunction(unifrnd(VarMin,VarMax,VarSize)));  % should be 2

%% MOEA/D Settings

MaxIt = 100;       % Maximum Number of Iterations
nPop  = 50;        % Population Size (Number of Sub-Problems)
nArchive = 50;     % External Archive Size

T = max(ceil(0.15*nPop),2);    % Number of Neighbors
T = min(max(T,2),15);

crossover_params.gamma  = 0.5;
crossover_params.VarMin = VarMin;
crossover_params.VarMax = VarMax;

%% Initialization

% Create Sub-problems
sp = CreateSubProblems(nObj,nPop,T);

% Empty Individual Template
empty_individual.Position = [];
empty_individual.Cost     = [];
empty_individual.g        = [];
empty_individual.IsDominated = [];

% Initialize Reference Point (ideal point)
z = zeros(nObj,1);

% Create Initial Population
pop = repmat(empty_individual,nPop,1);
for i=1:nPop
    % Sample respecting VarMin/VarMax
    pop(i).Position = VarMin + rand(VarSize).*(VarMax-VarMin);
    pop(i).Cost     = CostFunction(pop(i).Position);
    z = min(z,pop(i).Cost);
end

for i=1:nPop
    pop(i).g = DecomposedCost(pop(i),z,sp(i).lambda);
end

% Determine Domination Status
pop = DetermineDomination(pop);

% Initialize External Archive (Pareto Front Estimate)
EP = pop(~[pop.IsDominated]);

%% Main Loop
for it=1:MaxIt
    for i=1:nPop
        
        % Reproduction (Crossover)
        K = randsample(T,2);
        j1 = sp(i).Neighbors(K(1));
        j2 = sp(i).Neighbors(K(2));
        
        p1 = pop(j1);
        p2 = pop(j2);
        
        y = empty_individual;
        y.Position = Crossover(p1.Position,p2.Position,crossover_params);
        
        % Enforce Bounds
        y.Position = max(y.Position,VarMin);
        y.Position = min(y.Position,VarMax);
        
        y.Cost = CostFunction(y.Position);
        z = min(z,y.Cost);
        
        % Update Neighbors
        for j=sp(i).Neighbors
            y.g = DecomposedCost(y,z,sp(j).lambda);
            if y.g <= pop(j).g
                pop(j) = y;
            end
        end
    end
    
    % Update Domination
    pop = DetermineDomination(pop);
    ndpop = pop(~[pop.IsDominated]);
    
    EP = [EP
          ndpop]; %#ok
    EP = DetermineDomination(EP);
    EP = EP(~[EP.IsDominated]);
    
    % Control archive size
    if numel(EP) > nArchive
        Extra = numel(EP) - nArchive;
        ToBeDeleted = randsample(numel(EP),Extra);
        EP(ToBeDeleted) = [];
    end
    
    % Plot EP
    figure(1);
    PlotCosts(EP);
    xlabel('f1'); ylabel('f2');
    title(['MOEA/D: Iteration ' num2str(it)]);
    pause(0.01);
    
    % Display Iteration Info
    disp(['Iteration ' num2str(it) ...
          ': Number of Pareto Solutions = ' num2str(numel(EP))]);
end

%% Results


% Original Pareto
nPoints = 200;
ff1 = linspace(0, 1, nPoints);   % for TP1, x1 in [0,1]
ff2 = 1 - ff1.^2;                % TP1 Pareto front formula

% Build matrix of true Pareto front (2 x nPoints)
PFtrue = [ff1; ff2];


P_metrics = [EP.Cost];      % example: solver's Pareto points

%% Compute metrics
% -> Generational Distance (GD), a performance indicator in multi-objective optimization that measures
%    how far an algorithm's solutions are from the true Pareto front; 
% -> Spread (delta): measures how evenly the obtained Pareto solutions are distributed along the front

[GD, Delta] = ComputeMetrics(P_metrics, PFtrue);

%S metrics 2xN, 2x1
% P = [EP.Cost];   % 2 x N
% 
% % Worst per objective
% worst_vals = max(P, [], 2);
% 
% % Add margin (say 10%)
% ref = worst_vals + 0.1*abs(worst_vals);
% 
% S = SMetric2D(P, ref);
% 

disp('MOEA/D: Final Pareto front statistics:');
EPC = [EP.Cost];
for j=1:nObj
    disp(['Objective #' num2str(j) ':']);
    disp(['   Min = ' num2str(min(EPC(j,:)))]);
    disp(['   Max = ' num2str(max(EPC(j,:)))]);
    disp([' Range = ' num2str(max(EPC(j,:))-min(EPC(j,:)))]);
    disp(['  St.D. = ' num2str(std(EPC(j,:)))]);
    disp(['  Mean = ' num2str(mean(EPC(j,:)))]);
    disp(' ');
end
fprintf('GD = %.4f, Spread = %.4f\n', GD, Delta);
%disp(['S metric(hypervolume) = ' num2str(S)]);
disp(' ');
