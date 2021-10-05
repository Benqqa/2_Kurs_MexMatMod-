
Err=[]
Err_b=[]
Con=[]
K_o1=[]
K_o2=[]
% hilb(8)
n=10
for i_i=1:14
    d=0
while true
A=rand(n,n)
[Q,R]=qr(A)
zeros(1,n)
for i=1:n
D(i,i)=10^(14)*i^(-i_i)
end
% D=eye(n)*10^(-15)
d=1;

A=Q*D*Q'
A_ist=A
for k=0:(n-1)
    M=A(1:k,1:k);
    d=d*det(M)
end
if d~=0
    break
end
end
%  while A~=A.' 
% A =rand(10)
% end
Con=[Con,cond(A)]
root=ones(n,1)
b=A*root
b_ist=b

%root=[1;1;1;1;1;1;1;1;1;1]
A_1=A;
% A=10*randn(10)
n=length(b)
if det(A)~=0
A=[A,b]

%пр€мой ход метода √ауса
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
%     k=k+1
end
b=A(:,n+1)
A=A(:,1:n)
%обратный ход метода √ауса
x=zeros(n,1);
x(n)=b(n)/A(n,n);
for k=n-1:-1:1
    x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
end
x
Err=[Err,norm(x-root)/norm(x)]
% b_f=A*x
% Con
% Err_b=[Err_b,norm(b_f-b_ist)/norm(b_f)]
% K_o1=[K_o1,(norm(x-root)/norm(x))/(cond(A_ist)*norm(b_f-b_ist)/norm(b_f))]
% Err_A=norm(x*b_ist'-A_ist)/norm(x*b_ist')
% K_o2=[K_o2,(norm(x-root)/norm(root))/(cond(A_ist)*Err_A)]
end

% 
% %пара
% disp("число обусловленности A "+cond(A_1))
% B=[1,1;
%     1,1.0001]
% disp("число обусловленности B "+cond(B))
% norm([1;1]-[0;2])/norm([1;1])
% norm([2;2.0001]-[2;2.0002])/norm([2;2.0001])%числр обусленности матрицы cont(A)
% %C=rand(5)
% C=hilb(8)%нехороша€ матрица))))
% x_c=ones(8,1)
% b_c=C*x_c
% x_cc=C\b_c
% norm(x_cc-x_c)/norm(x_c)
% cond(C)
% 


figure
loglog(sort(Con),sort(Err),'c*')
% plot(Con,Err, 'c*')

% legend('Con','Err')
xlabel('Con')
ylabel('Err')

Err
Con

sort(Err);
sort(Con);


