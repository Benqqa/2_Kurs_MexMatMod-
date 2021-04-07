f= @(x)(2-x.^2).*cos(x)+2.*x.*sin(x)
syms x

a=-2
b=2
n=5

Xi_Ravn=GridRavn(a,b,n)
Xi_Cheb=GridCheb(a,b,n)
Poly_Spline_Cheb=Spline(a,b,Xi_Cheb,f,n)
Graphic_Spline(Poly_Spline_Cheb,Xi_Cheb,n,f,a,b)
Poly_Spline_Ravn=Spline(a,b,Xi_Ravn,f,n)
Graphic_Spline(Poly_Spline_Ravn,Xi_Ravn,n,f,a,b)

function [Grid]= GridRavn(a,b,n);
    Grid=[];
    h1=(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h1*i];
    end
end
function [Grid]=GridCheb(a,b,n);
    Grid=[];
    for i=0:n
        Grid=[Grid,0.5*(a+b)+0.5*(b-a)*cos((2*i+1)*pi/(2*(n+1)))];
    end
    Grid=sort(Grid)
end
function [M]=M_values(n,h,y,a,b)
    syms x
    func=(2-x.^2).*cos(x)+2.*x.*sin(x);
    df=diff(func);
    ddf=diff(df);
    M0=subs(ddf,x,a);
    Mn=subs(ddf,x,b);
    A=zeros(n+1);
    b=ones(n+1,1);
    for i=2:n
        A(i,i-1)=h(i)/(h(i)+h(i+1));
        A(i,i)=2;
        A(i,i+1)=h(i+1)/(h(i)+h(i+1));
    end
    for i=2:n
        b(i)=(6/(h(i)+h(i+1)))*((y(i+1)-y(i))/h(i+1)-(y(i)-y(i-1))/h(i));
    end
    b(1)=M0;
    b(n+1)=Mn;
    A(1,1)=1;
    A(1,2)=0;
    A(n+1,n+1)=1;
    A(n+1,n)=0;
    M=A\b;
end
function [g] = Spline(a,b,X,f,n);
    syms x;
    g=[];
    h=ones(n+1,1);
    y=ones(n+1,1);
    for i=2:n+1
        h(i)=X(i)-X(i-1);
    end
    for i=1:n+1
        y(i)=f(X(i));
    end
    M=M_values(n,h,y,a,b);
    for i=2:n+1
        gi=M(i-1)*(((X(i)-x)^3)/(6*h(i)))+M(i)*((x-X(i-1))^3)/(6*h(i))+(x-X(i-1))*(((y(i)-y(i-1))/h(i))-(h(i)/6)*(M(i)-M(i-1)))+(y(i-1)-M(i-1)*((h(i))^2)/6);
        g=[g,gi];
    end
end
function []=Graphic_Spline(Poly_mass,X,n,f,a,b);
    syms x;
    figure;
    hold on;
    plot(X,f(X),'o');
    plot(a:0.0001:b,f(a:0.0001:b),'b');
    for i=1:n
        X1=[X(i):0.001:X(i+1)];
        gi=Poly_mass(i);
        plot(X1,subs(gi, x, X1),'r');
    end
    title('График исходной функции и сплайна');
    legend('Узлы','Исходная функция','Cплайн');
    ylabel('Y');
    xlabel('X');
end