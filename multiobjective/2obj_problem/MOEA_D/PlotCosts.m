%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA124
% Project Title: Implementation of MOEA/D
% Muti-Objective Evolutionary Algorithm based on Decomposition
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function PlotCosts(EP)

    EPC=[EP.Cost];
    plot(EPC(1,:),EPC(2,:),'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 5,'DisplayName','non-dominating solutions');
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    grid on;

%plot the original effective pareto
%% original efficient front
hold on

ff1 = linspace(-1,1,200);   
ff2 = 1 - ff1.^2;            
plot(ff1, ff2, 'k-', 'LineWidth',1, 'DisplayName','original efficient front');

hold off;


end