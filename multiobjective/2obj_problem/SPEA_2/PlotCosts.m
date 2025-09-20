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
    plot(PFC(1,:),PFC(2,:),'ro');
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