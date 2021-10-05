function y = forwsub(A,b)
n = length(b); 
y = zeros(n,1); 
y(1) = b(1)/A(1,1);
for k = 2:n
    y(k)=(b(k) - A(k,1:k-1)*y(1:k-1))/A(k,k);
end