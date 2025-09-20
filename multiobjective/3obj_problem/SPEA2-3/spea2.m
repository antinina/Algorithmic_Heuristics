%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YOEA122
% Project Title: Strength Pareto Evolutionary Algorithm 2 (SPEA2)
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

% CostFunction=@(x) ZDT(x);
% 
% nVar=30;             % Number of Decision Variables
% 
% VarSize=[nVar 1];   % Decision Variables Matrix Size
% 
% VarMin=0;           % Decision Variables Lower Bound
% VarMax=1;           % Decision Variables Upper Bound
% 
% clc;
% clear;
% close all;

%% Problem Definition

nVar = 3;             % Number of Decision Variables
VarSize=[nVar 1];   % Decision Variables Matrix Size

% Variable bounds
VarMin = [0 -1 -1];
VarMax = [1  1  1];

% Cost function with parameters alpha=0.75, beta=10
CostFunction = @(x) three_obj_tp1(x,0.75,10);

%% SPEA2 Settings

MaxIt = 200;          % Maximum Iterations
nPop = 50;            % Population Size
nArchive = 50;        % Archive Size
K = round(sqrt(nPop+nArchive));  % KNN parameter

pCrossover = 0.7;
nCrossover = round(pCrossover*nPop/2)*2;
pMutation = 1-pCrossover;
nMutation = nPop - nCrossover;

crossover_params.gamma = 0.1;
crossover_params.VarMin = VarMin;
crossover_params.VarMax = VarMax;

mutation_params.h = 0.2;
mutation_params.VarMin = VarMin;
mutation_params.VarMax = VarMax;

%% Initialization

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.S = [];
empty_individual.R = [];
empty_individual.sigma = [];
empty_individual.sigmaK = [];
empty_individual.D = [];
empty_individual.F = [];

pop = repmat(empty_individual,nPop,1);
for i=1:nPop
    pop(i).Position = VarMin + (VarMax - VarMin) .* rand(1,3);

    pop(i).Cost = CostFunction(pop(i).Position);
end

archive = [];

%% Main Loop

for it=1:MaxIt
    
    % Combine population and archive
    Q = [pop; archive];
    nQ = numel(Q);
    
    % Reset strength
    for i=1:nQ
        Q(i).S = 0;
    end
    
    % Calculate strength
    dom = false(nQ,nQ);
    for i=1:nQ
        for j=i+1:nQ
            if Dominates(Q(i), Q(j))
                Q(i).S = Q(i).S + 1;
                dom(i,j) = true;
            elseif Dominates(Q(j), Q(i))
                Q(j).S = Q(j).S + 1;
                dom(j,i) = true;
            end
        end
    end
    
    % Raw fitness
    S = [Q.S];
    for i=1:nQ
        Q(i).R = sum(S(dom(:,i)));
    end
    
    % Density estimation
    Z = [Q.Cost]';
    SIGMA = pdist2(Z,Z,'seuclidean');
    SIGMA = sort(SIGMA);
    for i=1:nQ
        Q(i).sigma = SIGMA(:,i);
        Q(i).sigmaK = Q(i).sigma(K);
        Q(i).D = 1/(Q(i).sigmaK + 2);
        Q(i).F = Q(i).R + Q(i).D;
    end
    
    % Environmental selection
    nND = sum([Q.R]==0);
    if nND <= nArchive
        F = [Q.F];
        [~,SO] = sort(F);
        Q = Q(SO);
        archive = Q(1:min(nArchive,nQ));
    else
        archive = Q([Q.R]==0);
        SIGMA = SIGMA(:, [Q.R]==0);
        k = 2;
        while numel(archive) > nArchive
            while min(SIGMA(k,:))==max(SIGMA(k,:)) && k<size(SIGMA,1)
                k = k+1;
            end
            [~,j] = min(SIGMA(k,:));
            archive(j) = [];
            SIGMA(:,j) = [];
        end
    end
    
    PF=archive([archive.R]==0); % Approximate Pareto Front
    
    % Plot Pareto Front
    figure(1);
    PlotCosts(PF);
    pause(0.01);
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Number of PF members = ' num2str(numel(PF))]);
    
    if it>=MaxIt
        break;
    end
    
     % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for c = 1:nCrossover/2
        p1 = BinaryTournamentSelection(archive, [archive.F]);
        p2 = BinaryTournamentSelection(archive, [archive.F]);
        [popc(c,1).Position, popc(c,2).Position] = Crossover(p1.Position,p2.Position,crossover_params);
        
        % Keep within bounds
        popc(c,1).Position = max(popc(c,1).Position, VarMin);
        popc(c,1).Position = min(popc(c,1).Position, VarMax);
        popc(c,2).Position = max(popc(c,2).Position, VarMin);
        popc(c,2).Position = min(popc(c,2).Position, VarMax);
        
        popc(c,1).Cost = CostFunction(popc(c,1).Position);
        popc(c,2).Cost = CostFunction(popc(c,2).Position);
    end
    popc = popc(:);
    
    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for m=1:nMutation
        p = BinaryTournamentSelection(archive, [archive.F]);
        popm(m).Position = Mutate(p.Position, mutation_params);
        
        % Keep within bounds
        popm(m).Position = max(popm(m).Position, VarMin);
        popm(m).Position = min(popm(m).Position, VarMax);
        
        popm(m).Cost = CostFunction(popm(m).Position);
    end
    
    % New population
    pop = [popc; popm];
    
end

%% Results

% Original Pareto

P_solver = [PF.Cost];   % 3 x 50

% Generate true Pareto front (3D)
nPoints = 50;
[f1_grid, f2_grid] = meshgrid(linspace(0,1,nPoints));
f3_grid = 2 - f1_grid.^2 - f2_grid.^2;

PFtrue = [f1_grid(:)'; f2_grid(:)'; f3_grid(:)'];  % 3 x (50*50)

% Compute metrics
[GD, Delta] = ComputeMetrics(P_solver, PFtrue);


disp(' ');

PFC = [PF.Cost];
disp('SPEA-2:Final Pareto Front Statistics:');
for j=1:size(PFC,1)
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
