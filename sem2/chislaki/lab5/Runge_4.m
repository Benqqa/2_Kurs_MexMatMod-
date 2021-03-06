f_z= @(x,y,z)tan(x).*z-3*y+sin(x)
% f_z= @(x,y,z)2*x.*z+2*y-4*x
f_y= @(x,y,z) z
%y=@(x) sin(x)
Por=2
a=0
b=pi/2
% b=1.5
x0=a
z0=1
y0=0
eps=10^(-1)
[Err, Iter]=RungeToEps(eps,f_z,a,b,y0,z0) 
Y(1)=y0
Z(1)=z0
n=1
% [X,Y,Z]=runge(f_z,a,b,n,y0,z0)
% figure
% hold on
% grid on
% plot(X,Y)
% plot(X,sin(X))
%     legend('????? Y','???????')
function [Err, Iter] = RungeToEps(eps,f,a,b,y0,z0) 
    Eps=[];
    Err=[];
    n=2;
    [X0,Y0,Z0]=runge(f,a,b,n,y0,z0);
    i=1;
    while true
        i=i+1;
        [X,Y,Z]=runge(f,a,b,n^i,y0,z0);
        dY=abs(max(Y)-max(Y0))
        Err=[dY,Err]
        if(dY<eps)
            break
        end
%         if i>1
%             break
%         end
        Y0=Y;
        Z0=Z;
        X0=X;
    end
    Iter=i
    figure
    hold on
    grid on
    loglog(X,f(X,Y,Z))
    hold on
    loglog(X,sin(X))
    legend('????? Y','???????')
end
function [X,Y,Z]= runge(f,a,b,n,y0,z0)
    Y=[]
    Z=[]
    Y(1)=y0
    Z(1)=z0
    X=GridRavn(a,b,n);
    for i= 1:length(X)-1
        h=abs((a-b)/n);
       % h=0.1;
        k1=Z(i);
        g1=f(X(i),Y(i),Z(i))
        
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
function [Grid]= GridRavn(a,b,n);
    Grid=[];
    h1=(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h1*i];
    end
end