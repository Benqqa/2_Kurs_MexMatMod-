f=@(x,y,z) 2*x*z+2*y-4*x;
true_f=@(x)x+exp(x.^2);
a=0;
b=1;
n=20;
eps=10^(-1);
k=0
while true
    k=k+2
    i=2^k;
    if(k==1)
       x = create_gridRavn(a,b,i); 
       [Y] = Diff(f,a,b,i,x);
    end
    
    xx = create_gridRavn(a,b,2*i);
    [YY] = Diff(f,a,b,2*i,xx);
    Error = [Error,abs(YY(i)-Y(i))/3]
      if Error(k)<eps
          disp(k)
          disp(Error(k))
          break
      end
      x=xx
      Y=YY
end

% for k = 1:n
%     i=2^k
%     x = create_gridRavn(a,b,i)
%     xx = create_gridRavn(a,b,2*i)
%     [Y] = Diff(f,a,b,i,x)
%     [YY] = Diff(f,a,b,2*i,xx)
%     Error = [Error,abs(max(YY(k))-max(Y(k)))/3]
% end

figure
hold on
grid on
plot(x,Y,'-*b','Linewidth',1.2)
plot(x,true_f(x),'-*r','Linewidth',1.2)
legend('Численное решение м.Эйлера-Коши','Точное решение')
ylabel('Y');
xlabel('X');
title('График полученного решения y(x)')
% 
% figure
% hold on
% grid on
% plot(1:n,Error,'-*g','Linewidth',1.2)
% title('График зависимости погрешности численного решения от числа разбиений')
% ylabel('Погрешность по правилу Рунге');
% xlabel('Число разбиений n');

function [Y] = Diff(f,a,b,n,x)
y1=[];
z1=[];
%начальные условия
y=[1];
z=[1];
Y=[1];
Z=[1];
    for i=2:n+1
        y(i)=Y(i-1)+(b-a)/n*Z(i-1);
        y1=[y,y(i)];
        z(i)=Z(i-1)+(b-a)/n*f(x(i-1),Y(i-1),Z(i-1));
        z1=[z,z(i)];
        
        Y(i)=Y(i-1)+(b-a)/(2*n)*(Z(i-1)+z1(i));
        Y=[Y,Y(i)];
        
        Z(i)=Z(i-1)+(b-a)/(2*n)*(f(x(i-1),Y(i-1),Z(i-1))+f(x(i-1),y1(i),z1(i)));
        Z=[Z,Z(i)];
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
