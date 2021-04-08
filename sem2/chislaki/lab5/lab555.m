% f=@(x,y,z) 2*x*z+2*y-4*x;
f= @(x,y,z)tan(x).*z-3*y+sin(x)
true_f=@(x)x+exp(x.^2)
true_f= @(x) sin(x)
a=0
b=pi/2
n=5
Error =[];
for k = 1:n
    i=2^k
    x = create_gridRavn(a,b,i)
    xx = create_gridRavn(a,b,2*i)
    [Y] = Diff(f,a,b,i,x)
    [YY] = Diff(f,a,b,2*i,xx)
    Error = [Error,abs(max(YY)-max(Y))/3]
end

figure
hold on
grid on
plot(x,Y,'-*b','Linewidth',1.2)
plot(x,true_f(x),'-*r','Linewidth',1.2)
legend('Численное решение м.Эйлера-Коши','Точное решение')
ylabel('Y');
xlabel('X');
title('График полученного решения y(x)')

figure
hold on
grid on
plot(1:n,Error,'-*g','Linewidth',1.2)
title('График зависимости погрешности численного решения от числа разбиений')
ylabel('Погрешность по правилу Рунге');
xlabel('Число разбиений n');

function [Y] = Diff(f,a,b,n,x)
y1=[]
z1=[]
%начальные условия
y=[0]
z=[1]
Y=[0]
Z=[1]
    for i=2:n+1
        y(i)=y(i-1)+(b-a)/n*z(i-1);
        z(i)=z(i-1)+(b-a)/n*f(x(i-1),y(i-1),z(i-1));
        y1=[y,y(i)]
        z1=[z,z(i)]
    end
    for i=2:n+1
        Y(i)=y(i-1)+(b-a)/(2*n)*(z(i-1)+z1(i));
        Z(i)=z(i-1)+(b-a)/(2*n)*(f(x(i-1),y(i-1),z(i-1))+f(x(i-1),y1(i),z1(i)));
        Y=[Y,Y(i)]
        Z=[Z,Z(i)]
    end
    Y(n+1)=[];
end
function [x]= create_gridRavn(a,b,n)
    x=[0];
    h =(b-a)/(n);
    for i=1:n
       x=[x,a+h*i];
    end
end
