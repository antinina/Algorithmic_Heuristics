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

% CostFunction=@(x) MOP4(x);      % Cost Function
% 
% nVar=3;             % Number of Decision Variables
% 
% VarSize=[1 nVar];   % Size of Decision Variables Matrix
% 
% VarMin=-5;          % Lower Bound of Variables
% VarMax= 5;          % Upper Bound of Variables

CostFunction=@(x) TP1(x);      % Cost Function

nVar=3;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=-1;          % Lower Bound of Variables
VarMax= 1;          % Upper Bound of Variables


% Number of Objective Functions
nObj=numel(CostFunction(unifrnd(VarMin,VarMax,VarSize)));


%% NSGA-II Parameters

MaxIt=100;      % Maximum Number of Iterations

nPop=50;        % Population Size

pCrossover=0.7;                         % Crossover Percentage
nCrossover=2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

pMutation=0.4;                          % Mutation Percentage
nMutation=round(pMutation*nPop);        % Number of Mutants

mu=0.02;                    % Mutation Rate

sigma=0.1*(VarMax-VarMin);  % Mutation Step Size


%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];
%FCost
pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
   % pop(i).Position=unifrnd(VarMin,VarMax,VarSize);    
   %  pop(i).Cost=CostFunction(pop(i).Position);

    %first elements different bounds
    pop(i).Position = zeros(VarSize);
    pop(i).Position(1) = unifrnd(0, VarMax);
    pop(i).Position(2:end) = unifrnd(VarMin, VarMax, [1, 2]);
    pop(i).Cost=CostFunction(pop(i).Position);
end

% Non-Dominated Sorting
[pop, F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop, F]=SortPopulation(pop);


%% NSGA-II Main Loop

for it=1:MaxIt
    
    disp('Crossover');
    % Crossover
    popc=repmat(empty_individual,nCrossover/2,2);
    for k=1:nCrossover/2
        
        i1=randi([1 nPop]);
        p1=pop(i1);
        
        i2=randi([1 nPop]);
        p2=pop(i2);
        
        %needs to be bounded
       % [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position);
       [y1, y2] =  Crossover(p1.Position,p2.Position);
       
       lb = [0 -1 -1];
       ub = [1 1 1];
       popc(k,1).Position = min(max(y1, lb), ub);
       popc(k,2).Position = min(max(y2, lb), ub);
       
      disp(popc(k,1).Position)
      disp(popc(k,2).Position)
       
        popc(k,1).Cost=CostFunction(popc(k,1).Position);
        popc(k,2).Cost=CostFunction(popc(k,2).Position);
        
    end
    popc=popc(:);
    
    disp('Mutation');
    % Mutation
    popm=repmat(empty_individual,nMutation,1);
    for k=1:nMutation
        
        i=randi([1 nPop]);
        p=pop(i);
        
        %popm(k).Position=Mutate(p.Position,mu,sigma);
        y = Mutate(p.Position,mu,sigma);
        lb = [0 -1 -1];
        ub = [1 1 1];
    

        popm(k).Position = min(max(y, lb), ub);
       disp(popm(k).Position)
        popm(k).Cost=CostFunction(popm(k).Position);
        
    end
    
      % Merge
    pop=[pop
         popc
         popm]; %#ok
     
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    pop=SortPopulation(pop);
    
    % Truncate
    pop=pop(1:nPop);
    
    % Non-Dominated Sorting
    [pop, F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop, F]=SortPopulation(pop);
    
    % Store F1
    F1=pop(F{1});
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    
    % Plot F1 Costs
    figure(1);
    PlotCosts(F1);
    pause(0.01);
    if it == MaxIt
    hold on
    
    %% original efficient front
    
    ff1 = linspace(-1,1,200);   
    ff2 = 1 - ff1.^2;            
    plot(ff1, ff2, 'k-', 'LineWidth',1, 'DisplayName','original efficient front');
    
    hold off;



    end    

end

%% Results
% Original Pareto
nPoints = 200;
ff1 = linspace(0, 1, nPoints);   % for TP1, x1 in [0,1]
ff2 = 1 - ff1.^2;                % TP1 Pareto front formula

% Build matrix of true Pareto front (2 x nPoints)
PFtrue = [ff1; ff2];


P_metrics = [F1.Cost];      % example: solver's Pareto points

%% Compute metrics
% -> Generational Distance (GD), a performance indicator in multi-objective optimization that measures
%    how far an algorithm's solutions are from the true Pareto front; 
% -> Spread (delta): measures how evenly the obtained Pareto solutions are distributed along the front

[GD, Delta] = ComputeMetrics(P_metrics, PFtrue);

 
disp('Final Pareto Front Statistics:');

NSG = [F1.Cost];
for j = 1:size(NSG,1)
    disp(['Objective #' num2str(j) ':']);
    disp(['      Min = ' num2str(min(NSG(j,:)))]);
    disp(['      Max = ' num2str(max(NSG(j,:)))]);
    disp(['    Range = ' num2str(max(NSG(j,:))-min(NSG(j,:)))]);
    disp(['    St.D. = ' num2str(std(NSG(j,:)))]);
    disp(['     Mean = ' num2str(mean(NSG(j,:)))]);
    disp(' ');
end

fprintf('GD = %.4f, Spread = %.4f\n', GD, Delta);
disp(' ');

