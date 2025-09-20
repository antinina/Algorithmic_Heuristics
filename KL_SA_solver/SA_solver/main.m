%% Kod u kom se implementira minimizacija cutsize
% koristeci SA algoritam (anneal.m funkcija preuzeta sa Matlabovog sajta)

% Polazni problem, matrica susedstva gde su tezine grana realni brojevi
W_real = [
	0 0 0 0.5 0 0.5 0 0 0 0;
	0 0 0 0.25 0.25 0 0.25 0.25 0 0;
	0 0 0 0 0.5 0 0 0.5 0 0;
	0.5 0.25 0 0 0.25 1.0 0.75 0.25 0 0;
	0 0.25 0.5 0.25 0 0 0.58 1.08 0 0.33;
	0.5 0 0 1.0 0 0 0.5 0 1.0 0;
	0 0.25 0 0.75 0.58 0.5 0 0.58 0.5 0.83;
	0 0.25 0.5 0.25 1.08 0 0.58 0 0 1.33;
	0 0 0 0 0 1.0 0.5 0 0 0.5;
	0 0 0 0 0.33 0 0.83 1.33 0.5 0];

% Polazna matrica samo nenulte vrednosti se menjaju sa 1
W = double(W_real ~= 0)



n = size(W,1);

%Prema zadatku inicijalna particija je {1,3,5,6,7,9}

x0 = [1,0,1,0,1,0,1,0,1,0];


%% Graficka reprezentacija
G = graph(W);
A = []; B = [];
for i = 1:10    
    if x0(i) == 1
        A(end+1) = i;
    elseif x0(i) == 0
        B(end+1) = i;
    end    

end    

%% Assign colors based on partitions
nodeColors = zeros(1,n); 
nodeColors(A) = 1; % red
nodeColors(B) = 2; % blue
%%  Plot graph
%figure(1);
figure
p = plot(G, 'Layout', 'force');  % force-directed layout
p.NodeCData = nodeColors;        % color nodes
colormap([1 0 0; 0 0 1]);        % red = 1, blue = 2
p.MarkerSize = 7;                % node size
p.LineWidth = 1.5;               % edge thickness
p.NodeLabel = arrayfun(@num2str, 1:n, 'UniformOutput', false);
colorbar off;
title('Polazni problem');
grid on;

%% Simulated Annealing

% Loss function (cutsize)
loss = @(x) cut_loss(x, W);
% 
% Options
options = anneal();                % get defaults
options.Verbosity = 2;             % reports everything
options.InitTemp = 5;              % starting temperature
options.CoolSched = @(T) 0.99*T;   % slower cooling
options.Generator = @(x) partition_generator(x); % generise nasumicne particije cvorova i pokusava da dodje do resenja

% Run SA
[bestPart, bestCut] = anneal(loss, x0, options);
disp('__________________________________________')
disp('Najbolje resenje:'), disp(bestPart)
disp(['Cut size = ', num2str(bestCut)])

%Graficka reprezentacija
G = graph(W);
A = []; B = [];
for i = 1:10
    if bestPart(i) == 1
        A(end+1) = i;
    elseif bestPart(i) == 0
        B(end+1) = i;
    end    

end    

%% Assign colors based on partitions
nodeColors = zeros(1,n); 
nodeColors(A) = 1; % red
nodeColors(B) = 2; % blue
%%  Plot graph
%figure(2);
figure
p = plot(G, 'Layout', 'force');  % force-directed layout
p.NodeCData = nodeColors;        % color nodes
colormap([1 0 0; 0 0 1]);        % red = 1, blue = 2
p.MarkerSize = 7;                % node size
p.LineWidth = 1.5;               % edge thickness
p.NodeLabel = arrayfun(@num2str, 1:n, 'UniformOutput', false);
colorbar off;
title('Resenje nakon primene SA');
grid on;
