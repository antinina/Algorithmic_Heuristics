clc;
clear;
close all;

%% Problem Definition

CostFunction = @(x) TP1(x);  % Use your TP1 function

nVar = 3;                     % Number of Decision Variables
VarSize = [nVar 1];           % Decision Variables Matrix Size

%Granice  - za x1 su drugacije
VarMin = [0; -1*ones(nVar-1,1)];   % Lower bounds
VarMax = [1;  1*ones(nVar-1,1)];   % Upper bounds

%% SPEA2 Settings

MaxIt = 100;           % Maximum Number of Iterations
nPop = 50;             % Population Size
nArchive = 50;         % Archive Size
K = round(sqrt(nPop+nArchive));  % KNN Parameter

pCrossover = 0.5;
nCrossover = round(pCrossover*nPop/2)*2;
pMutation = 1 - pCrossover;
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

pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    pop(i).Position = VarMin + rand(VarSize).*(VarMax - VarMin);  % Respect bounds
    pop(i).Cost = CostFunction(pop(i).Position);
end

archive = [];

%% Main Loop

for it = 1:MaxIt
    
    Q = [pop; archive];
    nQ = numel(Q);
    
    dom = false(nQ, nQ);
    for i = 1:nQ
        Q(i).S = 0;
    end
    
    for i = 1:nQ
        for j = i+1:nQ
            if Dominates(Q(i), Q(j))
                Q(i).S = Q(i).S + 1;
                dom(i,j) = true;
            elseif Dominates(Q(j), Q(i))
                Q(j).S = Q(j).S + 1;
                dom(j,i) = true;
            end
        end
    end
    
    S = [Q.S];
    for i = 1:nQ
        Q(i).R = sum(S(dom(:,i)));
    end
    
    Z = [Q.Cost]';
    SIGMA = pdist2(Z,Z,'seuclidean');
    SIGMA = sort(SIGMA);
    for i = 1:nQ
        Q(i).sigma = SIGMA(:,i);
        Q(i).sigmaK = Q(i).sigma(K);
        Q(i).D = 1/(Q(i).sigmaK + 2);
        Q(i).F = Q(i).R + Q(i).D;
    end
    
    nND = sum([Q.R] == 0);
    if nND <= nArchive
        F = [Q.F];
        [~, SO] = sort(F);
        Q = Q(SO);
        archive = Q(1:min(nArchive, nQ));
    else
        SIGMA = SIGMA(:, [Q.R]==0);
        archive = Q([Q.R]==0);
        
        k = 2;
        while numel(archive) > nArchive
            while min(SIGMA(k,:)) == max(SIGMA(k,:)) && k < size(SIGMA,1)
                k = k + 1;
            end
            [~, j] = min(SIGMA(k,:));
            archive(j) = [];
            SIGMA(:,j) = [];
        end
    end
    
    PF = archive([archive.R]==0);  % Approximate Pareto Front
    
    % Plot Pareto Front
    figure(1);
    PlotCosts(PF);
    xlabel('f1'); ylabel('f2');
    title(['SPEA-2: Iteration ' num2str(it)]);
    pause(0.01);
    
    % Display Iteration Information
    disp(['SPEA-2:Iteration ' num2str(it) ': Number of PF members = ' num2str(numel(PF))]);
    
    % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for c = 1:nCrossover/2
        p1 = BinaryTournamentSelection(archive,[archive.F]);
        p2 = BinaryTournamentSelection(archive,[archive.F]);
        
        [popc(c,1).Position, popc(c,2).Position] = Crossover(p1.Position, p2.Position, crossover_params);
 
  %Ovo menjam jer imam granice za nezavisne promenljive
        
        popc(c,1).Position = max(popc(c,1).Position, VarMin);
        popc(c,1).Position = min(popc(c,1).Position, VarMax);
        popc(c,2).Position = max(popc(c,2).Position, VarMin);
        popc(c,2).Position = min(popc(c,2).Position, VarMax);
  
 %Ovo dalje ok      
        popc(c,1).Cost = CostFunction(popc(c,1).Position);
        popc(c,2).Cost = CostFunction(popc(c,2).Position);
    end
    popc = popc(:);
    
    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for m = 1:nMutation
        p = BinaryTournamentSelection(archive,[archive.F]);
        popm(m).Position = Mutate(p.Position, mutation_params);
        popm(m).Position = max(popm(m).Position, VarMin);
        popm(m).Position = min(popm(m).Position, VarMax);
        popm(m).Cost = CostFunction(popm(m).Position);
    end
    
    % Create New Population
    pop = [popc; popm];
    
end

%% Results
% Original Pareto
nPoints = 200;
ff1 = linspace(0, 1, nPoints);   % for TP1, x1 in [0,1]
ff2 = 1 - ff1.^2;                % TP1 Pareto front formula

% Build matrix of true Pareto front (2 x nPoints)
PFtrue = [ff1; ff2];


P_metrics = [PF.Cost];      % example: solver's Pareto points

%% Compute metrics
% -> Generational Distance (GD), a performance indicator in multi-objective optimization that measures
%    how far an algorithm's solutions are from the true Pareto front; 
% -> Spread (delta): measures how evenly the obtained Pareto solutions are distributed along the front

[GD, Delta] = ComputeMetrics(P_metrics, PFtrue);


disp('Final Pareto Front Statistics:');

PFC = [PF.Cost];
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

