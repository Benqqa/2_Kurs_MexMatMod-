n=10
Con=[]
Err=[]
Err2=[]
C=10^(-16)
 for i_i=1:16
%%%%встроенная функция
% C=C*10
% if C>=1
% C=C/10^(16)
% end
% A=sprandn(n,n,0.9,C)
%%%
%%% через диаганьную
A=randn(n,n);
% A(1,1)=1
% A(1,2)=1
% A(2,2)=2
% A(2,1)=1

[Q,R]=qr(A);
for i=1:n
D(i,i)=10^(16)*i^(-i_i);
end
A=Q*D*Q'
%%%
A
Con=[Con,cond(A)];
root=ones(n,1)
b=A*root
Q=zeros(n,n);
R=zeros(n,n);
if det(A)~=0
for j=1:n
    disp("шаг j="+j+" Q(:,j)=")
    Q(:,j)=A(:,j)
    
    for i=1:(j-1)
        disp("шаг j="+j+"шаг i="+i+" R(i,j)=")
        R(i,j)= Q(:,i)'*A(:,j)
        disp("шаг j="+j+"шаг i="+i+" Q(:,j)=")
        Q(:,j)=Q(:,j)- R(i,j)*Q(:,i)
        
    end
    disp("шаг j="+j+" R(j,j)=")
    R(j,j)=norm(Q(:,j),2)
    
    if R(j,j)==0
        disp(A(:,i)+"(ai линейно зависит от a1, . . . , ai?1")
        break
    end
    disp("шаг j="+j+" Q(:,j)=")
    Q(:,j)=Q(:,j)/R(j,j)
    
end
else
    disp("det=0")
end
Q
R
%обратный ход метода Гауса
x=zeros(n,1);
b=Q'*b
x(n)=b(n)/R(n,n);
for k=n-1:-1:1
    x(k)=(b(k)-R(k,k+1:n)*x(k+1:n))/R(k,k);
end
x
Err2=[Err2,norm(x-root)/norm(x)]
E=max(max(eye(n)-Q*Q'))

Err=[Err,max(max(eye(n)-Q*Q'))]
[Q1,R1]=qr(A);

 end
 figure
loglog(Con,Err,'c*')
% plot(Con,Err, 'c*')

% legend('Con','Err')
xlabel('Con')
ylabel('Err')
figure
loglog(Con,Err2,'c*')
% plot(Con,Err, 'c*')

% legend('Con','Err')
xlabel('Con')
ylabel('Err2')
