n=3
X=randn(n,n)
Y=zeros(n,n)
Y(:,1)=X(:,1)
for i=2:n
    i=i
    Y(:,i)=X(:,i)
    for j=1:(i-1)
        j=j
        ((X(:,i).*Y(:,j))./(Y(:,j).*Y(:,j))).*Y(:,j)
    Y(:,i)=Y(:,i)-((X(:,i).*Y(:,j))./(Y(:,j).*Y(:,j))).*Y(:,j)
    end
end
A=eye(n)
for i=2:n
    for j=1:(i-1)
        A(j,i)=(X(:,j) 
        for k=1:i
        A(j,i)=A(j,i)*X(:,k)
        end
    end
end