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
            for k=i:n+1
                t=A(i,k);
                A(i,k)=(c*A(i,k)+s*A(j,k));
                A(j,k)=(-s*t+c*A(j,k));
            end
        end
    end
    b=A(:,n+1)
    A=A(:,1:n)
    %обратный ход метода Гауса
    x=zeros(n,1);
    x(n)=b(n)/A(n,n);
    for k=n-1:-1:1
        x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
    end
    Err=[Err,norm(x-root)/norm(x)];
    x
    A=A_ist
    b=b_ist
    for k=0:(n-1)
            M=A(1:k,1:k);
            d=d*det(M)
        end
        if d~=0
            break
        end
    %СЛУЧАЙ ОПРЕДЕЛИТЕЛЯ=0
            if det(A)~=0
            %Объединяем в расширенную матрицу    
            A=[A,b]


            %Прямой ход метода Гаусса
                for i=1:(n)%Перебор строк
                    for j=(i+1):(n) %Перебор столбцов

            %Проверка на нулевой элемент
                        if A(i,i)==0
                            disp("ноль на диагонали")
                            max=0; %макс элемент
                            index_max=0; %индекс макс
                            for r=i+1:n %перебор строк
                                if abs(A(r,i))>max %если эл-т больше максимума
                                    max=abs(A(r,i)); %то он становится максимумом
                                    index_max=r; %забиваем его индекс
                                end
                            end
                            c=A(i,:);%берем i-тую строку и до конца 
                            A
                            A(i,:)=A(r,:);%меняем ее с r-той строкой,т.е максимальной
                            A(r,:)=c
                        end
            %Проверка закончена

                        m=A(j,i)/A(i,i);
                        for k=1:(n+1)%с вектором б, все элементы j-той строки
                            A(j,k)=A(j,k)-m*A(i,k); %берем этот элемент и вычитаем верхний
                        end
                    end
                end

            %СЛУЧАЙ ОПРЕДЕЛИТЕЛЯ!=0
            else
                disp("det!=0") 
            end
            %Разделяем расширенную матрицу
            b=A(:,n+1)
            A=A(:,1:n)

            %обратный ход метода Гаусса
            x=zeros(n,1);
            x(n)=b(n)/A(n,n);%найдем снизу справа последний
            for k=n-1:-1:1 %идем наверх для оставшихся
                x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k); %элемент с вектора переносим, умножая на предыдущий и на отношение
            end
            
   err=[err,(norm(x-root))/norm(x)]
end

figure
loglog(Con,Err,'r--o',Con,err,'g--o')
legend('Rotation', "Gauss")
xlabel('Con')
ylabel('Err')

