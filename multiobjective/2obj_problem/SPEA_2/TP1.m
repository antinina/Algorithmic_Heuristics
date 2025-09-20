function f = TP1(x)
    % Parameters
    a = 1;    
    b = 1;

    % Decision variables
    x1 = x(1); 
    h = 1 - x1^2;

    % For variables starting from index 2
    if length(x) > 1
        g = sum(10 + x(2:end).^2 - 10*cos(4*pi*x(2:end)));
    else
        g = 0;
    end

    %
    S = a/(0.2 + x1) + b*x1^2;    

    % Objectives
    f1 = x1;         % objective 1
    f2 = h + g*S;    % objective 2

    
    f = [f1
        f2];
end