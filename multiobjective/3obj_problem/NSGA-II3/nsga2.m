%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA120
% Project Title: Non-dominated Sorting Genetic Algorithm II (NSGA-II)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;

%% Problem Definition

nVar = 3;             % Number of Decision Variables
VarSize = [1 nVar];   % Size of Decision Variables Matrix

% Variable bounds
VarMin = [0 -1 -1];   % x1 ∈ [0,1], x2,x3 ∈ [-1,1]
VarMax = [1  1  1];

% Cost function with parameters alpha=0.75, beta=10
CostFunction = @(x) three_obj_tp1(x, 0.75, 10); 

% Number of Objectives
nObj = numel(CostFunction(unifrnd(VarMin,VarMax,VarSize)));

%% NSGA-II Parameters

MaxIt = 200;       % Maximum Number of Iterations
nPop = 50;         % Population Size

pCrossover = 0.7;                         
nCrossover = 2*round(pCrossover*nPop/2);  % Number of Offsprings

pMutation = 0.4;                          
nMutation = round(pMutation*nPop);        % Number of Mutants

mu = 0.02;                                % Mutation Rate
sigma = 0.1*(VarMax - VarMin);            % Mutation Step Size

%% Initialization

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.CrowdingDistance = [];

pop = repmat(empty_individual, nPop, 1);

for i = 1:nPop
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    pop(i).Cost = CostFunction(pop(i).Position);
end

% Non-Dominated Sorting
[pop, F] = NonDominatedSorting(pop);

% Crowding Distance
pop = CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F] = SortPopulation(pop);

%% NSGA-II Main Loop

for it = 1:MaxIt
    
    % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for k = 1:nCrossover/2
        i1 = randi([1 nPop]);
        p1 = pop(i1);

        i2 = randi([1 nPop]);
        p2 = pop(i2);

        [popc(k,1).Position, popc(k,2).Position] = Crossover(p1.Position, p2.Position);

        % Ensure bounds
        popc(k,1).Position = max(popc(k,1).Position, VarMin);
        popc(k,1).Position = min(popc(k,1).Position, VarMax);

        popc(k,2).Position = max(popc(k,2).Position, VarMin);
        popc(k,2).Position = min(popc(k,2).Position, VarMax);

        % Evaluate Cost
        popc(k,1).Cost = CostFunction(popc(k,1).Position);
        popc(k,2).Cost = CostFunction(popc(k,2).Position);
    end
    popc = popc(:);

    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for k = 1:nMutation
        i = randi([1 nPop]);
        p = pop(i);

        popm(k).Position = Mutate(p.Position, mu, sigma);

        % Ensure bounds
        popm(k).Position = max(popm(k).Position, VarMin);
        popm(k).Position = min(popm(k).Position, VarMax);

        % Evaluate Cost
        popm(k).Cost = CostFunction(popm(k).Position);
    end

    % Merge populations
    pop = [pop; popc; popm];

    % Non-Dominated Sorting
    [pop, F] = NonDominatedSorting(pop);

    % Crowding Distance
    pop = CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop, F] = SortPopulation(pop);

    % Truncate
    pop = pop(1:nPop);

    % Non-Dominated Sorting again
    [pop, F] = NonDominatedSorting(pop);

    % Crowding Distance
    pop = CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop, F] = SortPopulation(pop);

    % Store F1 (Pareto Front)
    F1 = pop(F{1});

    % Show Iteration Info
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);

    % Plot F1
    figure(1);
    PlotCosts(F1);
    pause(0.01);
end



%results
%% Results
% Original Pareto

P_solver = [F1.Cost];   % 3 x 50

% Generate true Pareto front (3D)
nPoints = 50;
[f1_grid, f2_grid] = meshgrid(linspace(0,1,nPoints));
f3_grid = 2 - f1_grid.^2 - f2_grid.^2;

PFtrue = [f1_grid(:)'; f2_grid(:)'; f3_grid(:)'];  % 3 x (50*50)

% Compute metrics
[GD, Delta] = ComputeMetrics(P_solver, PFtrue);

disp('Final Pareto Front Statistics:');

PFC = [F1.Cost];
for j = 1:size(PFC,1)
    disp(['Objective #' num2str(j) ':']);
    disp(['      Min = ' num2str(min(PFC(j,:)))]);
    disp(['      Max = ' num2str(max(PFC(j,:)))]);
    disp(['    Range = ' num2str(max(PFC(j,:))-min(PFC(j,:)))]);
    disp(['    St.D. = ' num2str(std(PFC(j,:)))]);
    disp(['     Mean = ' num2str(mean(PFC(j,:)))]);
    disp(' ');
end

fprintf('GD = %.4f, Spread = %.4f\n', GD, Delta);
disp(' ');

