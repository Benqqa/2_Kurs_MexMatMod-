fun= @(x) exp(-x.^2)
% fun= @(x) cos(x)
syms x;
a=-10
b=10
n=20

gridR=create_gridRavn(a,b,n)
% gridC=create_gridCheb(a,b,n)

% polySpline1=cubic_spline(a,b,gridC,fun,n)
% graphSpline(polySpline1,gridC,n,fun,a,b)

polySpline=cubic_spline(a,b,gridR,fun,n)
graphSpline(polySpline,gridR,n,fun,a,b)

function []=graphSpline(Poly_mass,grid,n,fun,a,b)
    syms x;
    figure
    hold on
    disp('график');
    plot(a:0.0001:b,fun(a:0.0001:b),'b')
    plot(grid,fun(grid),'*')
    for i =1:n
        X=[grid(i):0.001:grid(i+1)];
        pol=Poly_mass(i);
        plot(X,subs(pol, x, X))
    end
    legend('Табличная функция','Узлы','Естественный сплайн')
end
function [Poly] = cubic_spline(a,b,grid,fun,n)
    syms x;
    Poly=[];
    H=zeros(n+1,1);
    Y=zeros(n+1,1);
    for k=2:n+1
        H(k)=grid(k)-grid(k-1);
    end
    for k=1:n+1
        Y(k)=fun(grid(k));
    end
    M=find_M2(n,H,Y);
    for i=2:n+1
        s1=M(i-1)*(((grid(i)-x)^3)/(6*H(i)))
        s2=M(i)*((x-grid(i-1))^3)/(6*H(i))
        c1=((Y(i)-Y(i-1))/H(i))-(H(i)/6)*(M(i)-M(i-1))
        c2=Y(i-1)-M(i-1)*((H(i))^2)/6
        pol=s1+s2+(x-grid(i-1))*c1+c2
        Poly=[Poly,pol];
    end
end
function [M]=find_M2(n,H,Y)
    A=zeros(n+1);
    b=zeros(n+1,1);
    for i=2:n
            A(i,i-1)=H(i)/(H(i)+H(i+1));
            A(i,i)=2;
            A(i,i+1)=H(i+1)/(H(i)+H(i+1));
    end
    A(1,1)=1;
    A(n+1,n+1)=1;
    for i=2:n
        m1=6/(H(i)+H(i+1));
        s1=(Y(i+1)-Y(i))/H(i+1);
        s2=(Y(i)-Y(i-1))/H(i);
        b(i)=m1*(s1-s2);
    end
%     M=A\b
M=M_progonki(A,b);
end
function [Grid]= create_gridRavn(a,b,n)
    Grid=[];
    h=(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h*i];
    end
end

function Grid=create_gridCheb(a,b,n)
    Grid=[];
    for i=0:n
        Grid=[Grid,0.5*(a+b)+0.5*(b-a)*cos((2*i+1)*pi/(2*(n+1)))];
    end
    Grid=sort(Grid);
end
function X= M_progonki(A,b)
    M=A
    V=b
    A=diag(M,-1); % - поддиагональ матрицы коэффициентов
    B=diag(M);    % - главная диагональ матрици коэффициентов
    C=diag(M,1);  % - наддиагональ матрицы коэффициентов
    D=V;          % - вектор правой части системы
    n=length(B);
    Cz=zeros(1,n);
    Dz=zeros(1,n);
    Cz(1)=C(1)/B(1);
    for j=1:n-2
        Cz(j+1)=C(j+1)/(B(j+1)-A(j)*Cz(j));
    end
    Dz(1)=D(1)/B(1);
    for j=1:n-1
        Dz(j+1)=(D(j+1)-A(j)*Dz(j))/(B(j+1)-A(j)*Cz(j));
    end
    X(n)=Dz(n);
    for j=n-1:-1:1
        X(j)=Dz(j)-Cz(j)*X(j+1);
    end
end