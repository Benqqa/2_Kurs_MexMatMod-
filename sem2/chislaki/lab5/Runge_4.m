f= @(x,y,dy)tan(x)*dy-3*y+sin(x)
%y=@(x) sin(x)
Por=2
a=0
b=pi/2
x0=a
y0=0
dy0=1
Y=[]
dY=[]
Z=[]
Y(1)=y0
dY(1)=dy0
Z(1)=dy0
n=10
[Y,Z]=runge(f,a,b,n,Y,dY,Z)
X=GridRavn(a,b,n)
 
figure
hold on
grid on
plot(X,Y)
plot(X,sin(X))

function [Y,Z]= runge(f,a,b,n,Y,dY,Z)
    X=GridRavn(a,b,n)
    for i= 1:length(X)-1
        h=(a-b)/n;
        
        k1=Z(i);
        g1=f(X(i),Y(i),Z(i))
        
        k2=Z(i)+h*g1/2;
        g2=f(X(i)+h/2,Y(i)+h*k1/2,Z(i)+h*g1/2);
        
        k3=Z(i)+h*g2/2;
        g3=f(X(i)+h/2,Y(i)+h*k2/2,Z(i)+h*g2/2);
        
        k4=Z(i)+h*g3;
        g4=f(X(i)+h,Y(i)+h*k3,Z(i)+h*g3);
        
        dy=h/6*(k1+2*k2+2*k3+k4);
        dz=h/6*(g1+2*g2+2*g3+g4)
        Y(i+1)=Y(i)+dy;
        Z(i+1)=Z(i)+dz
    end
end
function [Grid]= GridRavn(a,b,n);
    Grid=[];
    h1=(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h1*i];
    end
end