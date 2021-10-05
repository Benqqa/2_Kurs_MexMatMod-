% A = [2,1,1;1,2,1;1,1,1]
% b=[4;4;3;4;4;3;4;4;3;4]
% b=[
%  47;
%  43;
%  43;
%  42;
% 121;
%  29;
%  49;
%  84;
%  23;
%   9]
Err=[]
Con=[]
% hilb(8)
n=10
C=10^(-16)
K=18
for i=1:K
    
%     C=-1*(3+randi(14))
C=C*10
if C>=1
C=C/10^(16)
end
% ran=randi(3)
% if ran==1
A=sprandn(n,n,1,C)
% else
%     A=rand(n,n)
% end
if i==K
    A=hilb(n)
end
% [Q,R,P]=qr(A)


%  while A~=A.' 
% A =rand(10)
% end
Con=[Con,cond(A)]
root=ones(n,1)
b=A*root

%root=[1;1;1;1;1;1;1;1;1;1]
A_1=A;
% A=10*randn(10)
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
%     k=k+1
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
Err=[Err,norm(x-root)/norm(x)]
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
% C=hilb(8)%нехорошая матрица))))
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
K

