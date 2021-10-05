n=10
Con=[]
Err=[]
err=[]
m=0
d=0
for i_i=1:16
    %%% генерация матрицы
    while true
        A=randn(n,n);
        [Q,R]=qr(A);
        for i=1:n
            D(i,i)=10^(16)*i^(-i_i);
        end
        A=Q*D*Q';
        a=A(1,1);
        b=A(2,1);
        c=a/sqrt(a^2+b^2);
        s=b/sqrt(a^2+b^2);
        %%% проверка условий применимости
        if (-s*a+c*b==0) && (c^2+s^2==1)
            break
        else
            disp("сгенерированая матрица не проходит по условию применимости, генерируем новую")
        end
    end
    A_ist=A
    Con=[Con,cond(A)];
    root=ones(n,1);
    %%%
    b=A*root% вектор b
    b_ist=b;
    A=[A,b];
    %%% метод вращений
    for i=1:n-1
        for j=i+1:n
            a=A(i,i);
            b=A(j,i);
            c=a/sqrt(a^2+b^2);
            s=b/sqrt(a^2+b^2);
            h=A(i,j)
            for k=i:n+1
                t=A(i,k);
                A(i,k)=(c*A(i,k)+s*A(j,k));
                A(j,k)=(-s*t+c*A(j,k));
%                 if i~=k
%                     if j~k
%                         A(i,j)=h
%                     end
%                 end
            end
        end
    end
    b=A(:,n+1)
    A=A(:,1:n)
    A
    %обратный ход метода Гауса
    x=zeros(n,1);
    x(n)=b(n)/A(n,n);
    for k=n-1:-1:1
        x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
    end
    Err=[Err,norm(x-root)/norm(x)];
    x
end

figure
loglog(sort(Con),sort(Err),'r--o')
legend('Rotation')
xlabel('Con')
ylabel('Err')

