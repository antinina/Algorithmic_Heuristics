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

function PlotCosts(PF)

    PFC=[PF.Cost];
      plot3(PFC(1,:),PFC(2,:),PFC(3,:),'r*','MarkerSize',8);
    xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
    
    title('SPEA-2: Non-dominated Solutions');
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