function cutVal = computeCut(W, A, B)
    %cutVal = sum(sum(W(A,B)));
    cutVal = nnz(W(A,B));  % nnz = number of nonzero elements
end
