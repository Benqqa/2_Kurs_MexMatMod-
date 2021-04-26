f= @(x,y,z)-1*z-x^2
f_ist = @(x) -x.^3/3+x.^2-2*x
%откуда

a=0;
b=1;
ya=f_ist(a)
yb=f_ist(b)

Iter=[];
Eps=[];

O1=-pi/3;
O2=pi/3;
% toch=5
% for i = 1:toch
%     
%     Error=[];
%     eps=10^(-i);
%     Eps=[Eps,eps];
%     k=0;
%     while true
%         k=k+1;
%         i=2^k;
%         if(k==1)
%             X0 = GridRavn(a,b,i);
%             [Y0] = m_shoot(f,a,b,i,X0,ya,yb,eps,O1,O2);
%         end
%         X = GridRavn(a,b,2*i);
%         [Y] = m_shoot(f,a,b,2*i,X,ya,yb,eps,O1,O2);
%         dY=abs(max(Y)-max(Y0));
%         Error = [Error,dY];
%           if dY<eps 
%               disp(k);
%               disp(Error(k));
%               break
%           end
%          X0=X;
%          Y0=Y;
%    end
%     Iter=[Iter,k]
% %     figure
% %     hold on
% %     grid on
% %     plot(X,Y,'Linewidth',3)
% %     plot([a:0.1:b],sin([a:0.1:b]),'Linewidth',1.5)
% %     legend('–унге','эталонное решение')
% 
%         
% end
% figure
% semilogx(Eps,Iter)
% title('зависимость числа итераций от заданной точности')
% figure
% semilogy(Error)
% title('зависимость ошибки от номера итерации')

n=2
X= GridRavn(a,b,n);
O1=-pi/3;
O2=pi/3;
[Y,O1,O2]=m_shoot(f,a,b,n,X,ya,yb,10^(-10),O1,O2)
% [Y,O1,O2]=m_shoot(f,a,b,n,X,ya,yb,0.00000000001,O1,O2)
figure
hold on
grid on

plot(X,Y,'Linewidth',3)
% plot([a:0.1:b],sin([a:0.1:b]),'Linewidth',1.5)
legend('эталонное решение','ѕристрелка','ѕристрелка n=4')

ErrAbs=[]
ErrRun=[]
for i=1:length(X)
   ErrAbs=[ErrAbs, abs(Y(i)-f_ist(X(i)))]
   ErrRun=[ErrRun, abs(Y(i)-f_ist(X(i)))/(2^2-1)]
end
figure
hold on
grid on

plot(ErrAbs,'Linewidth',3)
plot(ErrRun,'Linewidth',1.5)
legend('Abs','Runge')
function [Y,O1,O2] =m_shoot(f,a,b,n,X,ya,yb,eps,O1,O2)
    Y_n_1=runge(f,a,b,n,X,ya,tan(O1))
    Y_n_2=runge(f,a,b,n,X,ya,tan(O2))
    k=0;
    while true
        k=k+1
        disp(Y_n_1(length(Y_n_1)))
        disp(Y_n_2(length(Y_n_2)))
        d1=abs(Y_n_1(length(Y_n_1))-yb)
        d2=abs(Y_n_2(length(Y_n_2))-yb)
        if abs(O1-O2)<eps
            break
        end
        if d1<d2
            O2=(O2+O1)/2
            Y_n_2=runge(f,a,b,n,X,ya,tan(O2))
        else
            O1=(O1+O2)/2
            Y_n_1=runge(f,a,b,n,X,ya,tan(O1))
        end
    end
    Y=Y_n_1;
end
function [Grid]= GridRavn(a,b,n);
    Grid=[];
    h1=(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h1*i];
    end
end

function [Y]= runge(f,a,b,n,X,y0,z0)
    Y=[];
    Z=[];
    Y(1)=y0;
    Z(1)=z0;
    for i= 1:length(X)-1
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