n=5;
Err1=[];
Err2=[];
Err=[];
NormRaz=[];
ErrB1=[];
ErrB=[];
ErrB2=[];
p=0;
err=0;
m=0;
b0=0
m=m+1;
P=[];
for i_i=1:15
    %%conds
    A=condss(n,i_i)
    x0=ones(n,1);
    b=A*x0;
    n=length(b);
    x=Gauss1(A,b);
    normRaz=(norm(x-x0))/norm(x0);
    p=cond(A);
    P=[P,p] ;
    NormRaz=[NormRaz,normRaz];
end

%%2я зависимость
iter=100;
D1=eye(n);
for i=1:n
    D1(i,i)=10^(16)*i^(1);
end
D2=eye(n);
for i=1:n
    D2(i,i)=10^(16)*i^(5);
end
[OtnR1,OtnB1] = graph2(n,D1,iter);
[OtnR2,OtnB2] = graph2(n,D2,iter);
   
%%Graph 
figure
loglog(P,NormRaz,'k*') %строим график
hold on
grid on
title('График зависимости нормы разности решений от числа обусловленности')

figure
loglog(OtnB1,OtnR1 ) %строим график
legend('хорошо обусловленна')
title('График зависимости относительной ошибки в решении от относительного возмущения правой части при двух числах обусловленности')

figure
loglog(OtnB2,OtnR2 ) %строим график
legend('плохо обусловленна')


% loglog() %строим график

grid on
title('График зависимости относительной ошибки в решении от относительного возмущения правой части при двух числах обусловленности')
function [x] = Gauss1(A, b)
 n=length(b);

        %СЛУЧАЙ ОПРЕДЕЛИТЕЛЯ=0
            if det(A)~=0
            %Объединяем в расширенную матрицу    
            A=[A,b];


            %%Gauss
                for i=1:(n)%Перебор столбцов
                    for j=(i+1):(n) %Перебор строк

            %Проверка на нулевой элемент(выбор оптимального элемента)
                        if A(i,i)==0
%                             disp("ноль на диагонали")
                            max=0; %макс элемент
                            index_max=0; %индекс макс
                            for r=i+1:n %перебор элементов в строке
                                if abs(A(i,r))>max %если эл-т больше максимума
                                    max=abs(A(i,r)); %то он становится максимумом
                                    index_max=r; %запоминаем его индекс
                                end
                            end
                            c=A(:,i);%берем i-тую столбец
                            A;
                            A(:,i)=A(:,r);%меняем его с r-тым столбцом,т.е максимальный
                            A(:,r)=c;
                        end
            %Проверка закончена

                        m=A(j,i)/A(i,i);
                        for k=1:(n+1)%с вектором б, все элементы j-той строки
                            A(j,k)=A(j,k)-m*A(i,k); %берем этот элемент и вычитаем верхний
                        end
                    end
                end

            
            else
                disp("det!=0") 
            end
            %Разделяем расширенную матрицу
            b=A(:,n+1);
            A=A(:,1:n);

            %%Gauss_rev
            x=zeros(n,1);
            x(n)=b(n)/A(n,n);%найдем снизу справа последний
            for k=n-1:-1:1 %идем наверх для оставшихся
                x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k); %элемент с вектора переносим, умножая на предыдущий и на отношение
            end
end
function A= condss(n,i_i)
                    A=rand(n,n);
                    [Q,R]=qr(A);
                    D=eye(n);
                    for i=1:n
                        D(i,i)=10^(16)*i^(-i_i);
                    end
                    cond(A);
                    A=Q*D*Q';
end
function [OtnR,OtnB] = graph2(n,D,iter)
        OtnR=[];
        OtnB=[];
        A=rand(n,n);
        [Q,R]=qr(A);
        A=Q*D*Q';
        disp('acooooon')
        disp(cond(A))
        x=ones(n,1);
        b=A*x;
        x=Gauss1(A,b);
        for i=0:iter
            xi=ones(n,1);
            r=rand(n,1);
            dbi=r/norm(r,1)*(i/(25*iter)+1/100)*norm(b,1);
            bi=b+dbi;
            xi=Gauss1(A,bi);
            OtnR=[OtnR,norm(xi-x,1)/norm(xi,1)];
            OtnB=[OtnB,norm(bi-b,1)/norm(bi,1)];
        end
end