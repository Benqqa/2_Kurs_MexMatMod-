f= @(x,y,z)tan(x).*z-3*y+sin(x)
a=0;
b=pi/2;

z0=1
y0=0

Iter=[];
Eps=[];
toch=13
for i = 1:toch
    Error=[];
    eps=10^(-i);
    Eps=[Eps,eps];
    k=0;
    while true
        k=k+1;
        i=2^k;
        if(k==1)
            X0 = GridRavn(a,b,i);
            [Y0] = runge_Adams4(f,a,b,i,X0,y0,z0);
        end
        X = GridRavn(a,b,2*i);
        [Y] = runge_Adams4(f,a,b,2*i,X,y0,z0);
        dY=abs(max(Y)-max(Y0));
        disp(dY)
        Error = [Error,dY];
          if dY<eps 
              disp(k);
              disp(Error(k));
              break
          end
         X0=X;
         Y0=Y;

        
    end
    Iter=[Iter,k]
    figure
    hold on
    grid on
    plot(X,Y,'Linewidth',3)
    plot([a:0.1:b],sin([a:0.1:b]),'Linewidth',1.5)
    legend('–унге','эталонное решение')
    
end
        figure
    semilogx(Eps,Iter)
    title('зависимость числа итераций от заданной точности')
    figure
    semilogy(Error)
    title('зависимость ошибки от номера итерации')
function [Grid]= GridRavn(a,b,n);
    Grid=[];
    h1=(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h1*i];
    end
end
function [Y]=runge_Adams4(f,a,b,n,X,y0,z0)
    if(n>4)
        %разгон
        [Y,Z]= runge(f,a,b,n,X,y0,z0,3);
        %метод јдомсона4
        for i=4:length(X)-1
            
            Z(i+1)=Z(i)+((b-a)/n)*(-9/24*f(X(i-3),Y(i-3),Z(i-3))+37/24*f(X(i-2),Y(i-2),Z(i-2))-59/24*f(X(i-1),Y(i-1),Z(i-1))+55/24*f(X(i),Y(i),Z(i)));
            Y(i+1)=Y(i)+((b-a)/n)*(-9/24*Z(i-3)+37/24*Z(i-2)-59/24*Z(i-1)+55/24*Z(i));
        end
    else
       [Y,Z]= runge(f,a,b,n,X,y0,z0,length(X)-1);
    end
    
end
function [Y,Z]= runge(f,a,b,n,X,y0,z0,ogr)
    Y=[];
    Z=[];
    Y(1)=y0;
    Z(1)=z0;
    for i= 1:ogr
        h=abs((b-a)/n);
        k1=Z(i);
        g1=f(X(i),Y(i),Z(i));
        k2=Z(i)+h*g1/2;
        g2=f(X(i)+h/2,Y(i)+h*k1/2,Z(i)+h*g1/2);      
        k3=Z(i)+h*g2/2;
        g3=f(X(i)+h/2,Y(i)+h*k2/2,Z(i)+h*g2/2);        
        k4=Z(i)+h*g3;
        g4=f(X(i)+h,Y(i)+h*k3,Z(i)+h*g3);  
        dy=h/6*(k1+2*k2+2*k3+k4);
        dz=h/6*(g1+2*g2+2*g3+g4);
        Y(i+1)=Y(i)+dy;
        Z(i+1)=Z(i)+dz;
    end
end