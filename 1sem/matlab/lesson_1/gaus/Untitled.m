n=10
C=10^(-16)
K=18
if C>=1
C=C/10^(16)
end
A=rand(n,n)
if i==K
    A=hilb(n)
end
root=ones(n,1)
b=A*root
A_1=A;
n=length(b)
if det(A)~=0
  %проверить на 0 на дивагонали до самого применения!!!! чкать лекуию
A=[A,b]

%прямой ход метода Гауса
    for i=1:(n)
        for j=(i+1):(n)
            if A(i,i)==0
                disp("0 на диагонали")
%                 k=k+1
                break
            end
            m=A(j,i)/A(i,i);
            for k=1:(n+1)
                A(j,k)=A(j,k)-m*A(i,k);
            end
        end
    end

else
    disp("det=0")
end
b=A(:,n+1)
A=A(:,1:n)
%обратный ход метода Гауса
x=zeros(n,1);
x(n)=b(n)/A(n,n);
for k=n-1:-1:1
    x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
end
x
