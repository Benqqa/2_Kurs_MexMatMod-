n=10
A = randn(n,n);
[n,m]=size(A);
if n~=m
    error('The given matrix is not square!');
end
[B{1:n}]=deal(A);
p=ones(1,n+1);
for k=2:n
    p(k)=-trace(B{k-1})/(k-1);
    B{k}=A*B{k-1}+A*p(k)*eye(n));
end
p(n+1)=-trace(B{n})/n;
poly(A)