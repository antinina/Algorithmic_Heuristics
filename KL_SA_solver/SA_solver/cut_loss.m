function c = cut_loss(x, A)
    n = length(x);

    % Particija ne moze biti svi cvorovi ili nijedan
    if all(x==0) || all(x==1)
        c = inf;
        return
    end

    % x binarni vektor, A matrica susedstva
    
    c = 0;
    for i = 1:n
        for j = i+1:n
            if A(i,j) == 1 && x(i) ~= x(j)
                c = c + 1; % edge is cut
            end
        end
    end
    
    fprintf(1,'Curent cut size = %d\n',c);    
    disp('Partition :')
    disp(x)

end
