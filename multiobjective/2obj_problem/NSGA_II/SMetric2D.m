function S = SMetric2D(P, ref)
%   P   = 2 x N matrix of non-dominated solutions (each column is a point)
%   ref = 2 x 1 reference point (worse than all solutions)
%
%   S  = scalar value of the hypervolume

    % Ensure input is 2 x N
    if size(P,1) ~= 2
        error('Only works for 2D objectives (2 x N matrix expected)');
    end
    
    % Sort solutions by f1 (first objective)
    [~, idx] = sort(P(1,:));
    P = P(:, idx);
    
    % Initialize
    S = 0;
    prev_f1 = ref(1);
    N = size(P,2);
    for i = N:-1:1
        f1 = P(1,i);
        f2 = P(2,i);
        
        width  = prev_f1 - f1;
        height = ref(2) - f2;
        
        S = S + width * height;
        
        prev_f1 = f1;
    end
end
