function x = backsub(A,b)
    n = length(b); 
    x = zeros(n,1); 
    x(n) = b(n)/A(n,n);
    for k = n-1:-1:1
        x(k) = (b(k) - A(k,k+1:n)*x(k+1:n))/A(k,k);
end