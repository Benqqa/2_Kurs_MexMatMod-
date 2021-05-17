fun= @(x) sqrt(sin(x.^2))
a1=-3
b1=3
a2=4
b2=5
n=20
ErrR_L1=[]
ErrR_L2=[]
ErrC_L1=[]
ErrC_L2=[]

ErrR_S1=[]
ErrR_S2=[]
ErrC_S1=[]
ErrC_S2=[]
%  figure
for i=1:n
    [errR_L1,errC_L1,errR_S1,errC_S1]=interpol(a1,b1,fun,i)
    [errR_L2,errC_L2,errR_S2,errC_S2]=interpol(a2,b2,fun,i)
    ErrR_L1=[ErrR_L1,errR_L1]
    ErrR_L2=[ErrR_L2,errR_L2]
    ErrC_L1=[ErrC_L1,errC_L1]
    ErrC_L2=[ErrC_L2,errC_L2]
    
    ErrR_S1=[ErrR_S1,errR_S1]
    ErrR_S2=[ErrR_S2,errR_S2]
    ErrC_S1=[ErrC_S1,errC_S1]
    ErrC_S2=[ErrC_S2,errC_S2]
end
figure
hold on
plot(2:n+1,ErrR_L1,2:n+1,ErrC_L1)
plot(2:n+1,ErrR_S1,2:n+1,ErrC_S1)
legend('равномерная сетка L', 'чебышевская сетка L','равномерная сетка S', 'чебышевская сетка S')
title('max отклонение интерполяции от числа узлов на [-3,3]')
figure
hold on
plot(2:n+1,ErrR_L2,2:n+1,ErrC_L2)
plot(2:n+1,ErrR_S2,2:n+1,ErrC_S2)
legend('равномерная сетка L', 'чебышевская сетка L','равномерная сетка S', 'чебышевская сетка S')
title('max отклонение интерполяции от числа узлов на [4,5]')

function [ErrR_L,ErrC_L,ErrR_S,ErrC_S]=interpol(a,b,fun,n)
gridR=create_gridRavn(a,b,n)
gridC=create_gridCheb(a,b,n)

polyR=polyL(fun,gridR,n)
polyC=polyL(fun,gridC,n)

polySplineC=cubic_spline(a,b,gridC,fun,n)
polySplineR=cubic_spline(a,b,gridR,fun,n)

syms z x;
ErrR_L=max(abs(subs(polyR, z, a:0.01:b)-fun(a:0.01:b)))
ErrC_L=max(abs(subs(polyC, z, a:0.01:b)-fun(a:0.01:b)))
 mR=[]
 mC=[]
for i=1:n
    mR=[mR,max(abs(subs(polySplineR(i), x, gridR(i):0.01:gridR(i+1))-fun(gridR(i):0.01:gridR(i+1))))]
    mC=[mC,max(abs(subs(polySplineC(i), x, gridC(i):0.01:gridC(i+1))-fun(gridC(i):0.01:gridC(i+1))))]
end
ErrR_S=max(mR)
ErrC_S=max(mC)

graphSpline(polySplineC,gridC,n,fun,a,b,1)
graphSpline(polySplineR,gridR,n,fun,a,b,0)
graphLagrang(fun,a,b,polyC,gridC)
end
function g= graphLagrang(fun,a,b,poly,grid1)
    syms z;
    figure
    hold on
    grid on
    plot(grid1,fun(grid1),'*')
    plot(a:0.01:b, subs(poly, z, a:0.01:b))
    plot(a:0.0001:b,fun(a:0.0001:b))

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
    Grid=sort(Grid)
end
function [pol] = polyL(fun,grid,n)
    syms z;
    syms m;
    pol=0;
    m=1;
    for i=1:n+1
         m=1;
        for j=1:n+1
            if i~=j
                 m=m*(z-grid(j))/(grid(i)-grid(j));
            end
        end
        pol=pol+fun(grid(i))*m;
    end
end
function []=graphSpline(Poly_mass,grid,n,fun,a,b,index)
    syms x;
    figure
    hold on
    disp('график');
    plot(a:0.0001:b,fun(a:0.0001:b),'b')
    plot(grid,fun(grid),'*')
    for i =1:n
        X=[grid(i):0.001:grid(i+1)];
        pol=Poly_mass(i);
        plot(X,subs(pol, x, X),'g')
    end
    legend('Табличная функция','Узлы','Естественный сплайн')
    if index==1
        title('Чебышевская сетка')
    end
    if index==0
        title('Равномерная сетка')
    end
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
        s1=M(i-1)*(((grid(i)-x)^3)/(6*H(i)));
        s2=M(i)*((x-grid(i-1))^3)/(6*H(i));
        c1=((Y(i)-Y(i-1))/H(i))-(H(i)/6)*(M(i)-M(i-1));
        c2=Y(i-1)-M(i-1)*((H(i))^2)/6;
        pol=s1+s2+(x-grid(i-1))*c1+c2;
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
    A(1,2)=0;
    A(n+1,n+1)=1;
    A(n+1,n)=0;
    for i=2:n
        m1=6/(H(i)+H(i+1));
        s1=(Y(i+1)-Y(i))/H(i+1);
        s2=(Y(i)-Y(i-1))/H(i);
        b(i)=m1*(s1-s2);
    end
%     M=A\b
    M=M_progonki(A,b)
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