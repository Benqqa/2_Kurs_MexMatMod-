f= @(x,y,z)tan(x).*z-3*y+sin(x)
f_ist = @(x)sin(x)
p=@(x) -2*x;
q=@(x) -2;
f_i=@(x) -4*x;

a=0;
b=1;
ya=f_ist(a);
yb=f_ist(b);

n=7;
% [X,h]= GridRavn(a,b,n)
Iter=[];
Eps=[];
toch=5
for i = 1:toch
    
    Error=[];
    eps=10^(-i);
    Eps=[Eps,eps];
    k=0;
    while true
        k=k+1;
        i=2^k;
        if(k==1)
            [X0,h0] = GridRavn(a,b,i);
            [Y0] = Kon_Raz(f,p,q,f_i,a,b,i,X0,h0,ya,yb);
        end
        [X,h] = GridRavn(a,b,2*i);
        [Y] = Kon_Raz(f,p,q,f_i,a,b,2*i,X,h,ya,yb);
%         dY=abs(max(Y)-max(Y0));
        dY=abs(Y(length(X)-2)-Y0(length(X0)-1));
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
    plot([a:0.1:b],f_ist([a:0.1:b]),'Linewidth',1.5)
    legend('sss','эталонное решение')

        
end
figure
semilogx(Eps,Iter)
title('зависимость числа итераций от заданной точности')
figure
semilogy(Error)
title('зависимость ошибки от номера итерации')
n=100
[X,h] = GridRavn(a,b,n);
[Y] = Kon_Raz(f,p,q,f_i,a,b,n,X,h,ya,yb);
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
function [Y] = Kon_Raz(f,p,q,f_i,a,b,n,X,h,ya,yb)

P=zeros(length(X),1);
Q=zeros(length(X),1);
F=zeros(length(X),1);
for i=1:length(X)
P(i)=p(X(i));
Q(i)=q(X(i));
F(i)=f_i(X(i))*h^2;
% F(i)=f_i(X(i));
end
P;
Q;
F;
F(1)=ya;
% F(1)=1
F(length(X))=yb;
% F(length(X))=1
W=zeros(length(X),length(X));
for i=1:length(X)
        W(i,i)=Q(i)*h^2-2;
%          W(i,i)=1
     if(i<length(X))
        W(i,i+1)=1+P(i)*h/2;
%         W(i,i+1)=2
     end
     if(i>1)
         W(i,i-1)=1-P(i)*h/2;
%          W(i,i-1)=-2
     end
end
% W(1,1)=ya;
W(1,1)=1;
W(1,1+1)=0;
% W(length(X),length(X))=yb;
W(length(X),length(X))=1;
W(length(X),length(X)-1)=0;
length(X);

M=W;
V=F;
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
Y(n)=Dz(n);
for j=n-1:-1:1
    Y(j)=Dz(j)-Cz(j)*Y(j+1);
end
Y   ;
if det(W) == 0
    disp('det=0')
end


end

function [Grid,h1]= GridRavn(a,b,n);
    Grid=[];
    h1=abs(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h1*i];
    end
end