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

function PlotCosts(pop)

    Costs=[pop.Cost];
    
    plot3(Costs(1,:),Costs(2,:),Costs(3,:),'r*','MarkerSize',8);
    xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
    
    title('NSGA-II: Non-dominated Solutions (F_{1})');
    grid on;box on;

    
  
hold on;    
    %% ovo za pareto koji je vec poznat
    % Define a grid over f1 and f2
f1_vals = linspace(0,1,50);
f2_vals = linspace(-1,1,50);
[F1_grid, F2_grid] = meshgrid(f1_vals, f2_vals);

% Compute f3 values
F3_grid = 2 - F1_grid.^2 - F2_grid.^2;

% Plot the 3D surface

surf(F1_grid, F2_grid, F3_grid,'DisplayName','original efficient front');
xlabel('f1'); ylabel('f2'); zlabel('f3');


hold off;

end