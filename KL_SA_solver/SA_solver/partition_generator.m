function newx = partition_generator(x)
    % Ensure equal partition sizes
    n = length(x);
    part1 = find(x==0);
    part2 = find(x==1);

    % Pick random nodes and swap
    i = part1(randi(length(part1)));
    j = part2(randi(length(part2)));

    newx = x;
    newx(i) = 1;
    newx(j) = 0;
end
