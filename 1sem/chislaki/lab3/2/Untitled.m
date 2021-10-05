n=10
Con=[]
Err_G=[]
Err_g=[]
Err2_g=[]
Err_g2=[]
Err2_g2=[]
C=10^(-16)
 for i_i=1 :16
%%%%встроенная функция
% C=C*10
% if C>=1
% C=C/10^(16)
% end
% A=sprandn(n,n,0.9,C)
%%%
%%% через диаганьную
% d=0
% while true
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
A_ist=A
Con=[Con,cond(A)];
root=ones(n,1)
b=A*root
b_ist=b
% for k=0:(n-1)
%     M=A(1:k,1:k);
%     d=d*det(M)
% end
% if d~=0
%     break
% end
% end

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
% %обратный ход метода Гауса
% x=zeros(n,1);
% b=Q'*b
% x(n)=b(n)/R(n,n);
% for k=n-1:-1:1
%     x(k)=(b(k)-R(k,k+1:n)*x(k+1:n))/R(k,k);
% end
x=inv(R)*Q'*b
x
Err_g=[Err_g,norm(x-root)/norm(root)]
E=max(max(eye(n)-Q*Q'))
Err2_g=[Err2_g,max(max(eye(n)-Q*Q'))]

% Err_g=[Err_g,max(max(eye(n)-Q*Q'))]
[Q1,R1]=qr(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  while A~=A.' 
% A =rand(10)
% end
% Con=[Con,cond(A)]
root=ones(n,1)
A=A_ist
b=A*root
b_ist=b

%root=[1;1;1;1;1;1;1;1;1;1]
A_1=A;
% A=10*randn(10)
n=length(b)
if det(A)~=0
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
Err_G=[Err_G,norm(x-root)/norm(root)]

%%%%%%%%%%%%%%%%%%%%
root=ones(n,1)
A=A_ist
b=A*root
b_ist=b
Q=zeros(n,n);
V=zeros(n,n);
R=zeros(n,n);
if det(A)~=0
for i = 1:n
     V(:,i)=A(:,i);
end
for i = 1:n
    R(i,i)=norm(V(:,i),2);
    Q(:,i)=V(:,i)./R(i,i);
    for j =(i+1):n
        R(i,j)=Q(:,i)'*V(:,j);
        V(:,j)=V(:,j)-R(i,j).*Q(:,i);
    end
 end
else
    disp("det=0")
end
Q
R
%обратный ход метода Гауса
x=zeros(n,1);
% b=Q'*b
% x(n)=b(n)/R(n,n);
% for k=n-1:-1:1
%     x(k)=(b(k)-R(k,k+1:n)*x(k+1:n))/R(k,k);
% end

x=inv(R)*Q'*b
Err_g2=[Err_g2,norm(x-root)/norm(root)]
Err2_g2=[Err2_g2,max(max(eye(n)-Q*Q'))]
 end
 
 
 
 
 
 
 figure
loglog(Con,Err_g,'c--o',Con,Err_G,'g--o',Con,Err_g2,'r--o')
legend("Grama-Schmidta","Gauss","Mod-Grama-Schmidta")
% plot(Con,Err, 'c*')

% legend('Con','Err')
xlabel('Con')
ylabel('Err')
figure
loglog(Con,Err2_g,'c*',Con,Err2_g2,'g*')
% plot(Con,Err, 'c*')
legend("Grama-Schmidta","Mod-Grama-Schmidta")
% legend('Con','Err')
xlabel('Con')
ylabel('Err2')
