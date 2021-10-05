function X = res(A)
    n = size(A, 1);
    X(n) = 0;
    for i = n:-1:1
        for k = i:n
            A(i, n + 1) = A(i, n + 1) - A(i, k) * X(k);
        end
        X(i) = A(i, n + 1) / A(i, i);
    end
end