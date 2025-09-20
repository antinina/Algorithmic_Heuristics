% Test
% W = [0 1 1 0 0 0;
%      1 0 1 1 0 0;
%      1 1 0 1 0 0;
%      0 1 1 0 1 1;
%      0 0 0 1 0 1;
%      0 0 0 1 1 0];

% Matrica susedstva kola iz primera
W = [
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

n = size(W,1);

% Inicijalne particije {acegi , bdfghj}
A = [1 3 5 7 9];
B = [2 4 6 8 10];

cutFirst = computeCut(W, A, B);

%prikazi pocetni problem


%% Create graph
G = graph(W);

%% Assign colors based on partitions
nodeColors = zeros(1,n); 
nodeColors(A) = 1; % red
nodeColors(B) = 2; % blue
%%  Plot graph
figure(1);
p = plot(G, 'Layout', 'force');  % force-directed layout
p.NodeCData = nodeColors;        % color nodes
colormap([1 0 0; 0 0 1]);        % red = 1, blue = 2
p.MarkerSize = 7;                % node size
p.LineWidth = 1.5;               % edge thickness
p.NodeLabel = arrayfun(@num2str, 1:n, 'UniformOutput', false);
colorbar off;
title('Pocetni problem');
grid on;

%-------primena algoritma za trazenje min cut-size-a

[A,B,cutVal] = kernighan_lin(W, A, B);

disp(['' ...
    '*****************REZULTAT**********************']);
disp('Prva particija:'), disp(A)
disp('Druga particija:'), disp(B)
disp(['Minimalan cut size = ', num2str(cutVal)])
disp(['Vrednost cut size-a pocetnog problema = ', num2str(cutFirst)])
disp('_______________________________________________');
%% Create graph
G = graph(W);

%% Assign colors based on partitions
nodeColors = zeros(1,n); 
nodeColors(A) = 1; % red
nodeColors(B) = 2; % blue

%%  Plot graph
figure;
p = plot(G, 'Layout', 'force');  % force-directed layout
p.NodeCData = nodeColors;        % color nodes
colormap([1 0 0; 0 0 1]);        % red = 1, blue = 2
p.MarkerSize = 7;                % node size
p.LineWidth = 1.5;               % edge thickness
p.NodeLabel = arrayfun(@num2str, 1:n, 'UniformOutput', false);
colorbar off;
title('Graf kola nakon primene KL algoritma');
grid on;
