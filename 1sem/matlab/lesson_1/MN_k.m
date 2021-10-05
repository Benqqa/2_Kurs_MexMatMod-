syms x;
%входные данне
% a=-2.73;
% b= -1.59;
a=-1;
b= 1;

 fun= x.^2;
%  a=-2;
% b=4;
%  fun=3*exp(x)+5*x+2;
%  a=-0.2
%  b=0.2
%  fun=x+x.^2
% fun =(x-0.7).^(2);
m=1;
%%%%%%%%%%


fun_d=diff(fun);
fun_d2=diff(fun_d);


f=inline(fun);
f_d=inline(fun_d);
f_d2=inline(fun_d2);

stap=0.00001;
X=[a:stap:b];
Y=f(X);
Y_d=f_d(X);
Y_d_2=f_d2(X);
Err=[];
M2=max(abs(Y_d_2))
m1=min(abs(Y_d))
eps=10^(-3);
root=-0.69842
Eps=[]
k_eps=[];
S_iter=[];
% while eps>10^(-15)
% eps=eps/10
%  Eps=[Eps,eps]
%  if f(a)*f(b)<0 

disp("корень есть на ["+a+","+b+"]")
%     if ((Y_d>0) | (Y_d*-1>0))
        disp("первая производная знакопостяонна на ["+a+","+b+"]")
%         if ((Y_d_2>0) | (Y_d_2*-1>0))
            disp("вторая производная знакопостяонна на ["+a+","+b+"]")
            x0=b;
            while f(x0)*f_d2(x0)<=0
                x0=x0-stap;
            end
            disp("x0 = "+x0)
            xi=x0;
            xi_1=0;
            k=0;
%             while 0.5*(M2/m1)*abs(xi-xi_1)>=eps
            while abs(xi-xi_1)>=eps
                Err=[Err,abs(xi-root)]
                %Err=[Err,abs(0.5*abs(xi-xi_1)-root)]
                k=k+1;
                xi_1=xi;
                disp("приближение = "+xi)
                
                xi=xi-m*(f(xi)/f_d(xi));
                k_eps=[k_eps,k];
                %disp("корень не найти ")
            end
           S_iter=[S_iter,k]
            
            disp("корень = "+xi)
            koren=xi;
%         end
%     end
% else
%     disp("корень не найти ")
% end
% end
% % % % % % % % % 
figure
 hold on
 grid on
 title("fun= x.^2;")
%  title(" y=3*exp(x)+5*x+2;")
 xlabel('число итераций')
 ylabel('ошибка')
 plot(k_eps,Err)
% % % % % % % % % 
% figure
%  hold on
%  grid on
%  title(" y=x.^4 - x.^3 - 2*x.^2 + 3*x-3;")
%  % title(" y=3*exp(x)+5*x+2;")
%  xlabel('Ox')
%  ylabel('Oy')
% plot(X,f(X))
% plot(koren,0,'c*')
% % % % % % % % % 
% figure
%  hold on
%  grid on
% %  title(" y=x.^4 - x.^3 - 2*x.^2 + 3*x-3;")
%  title(" y=3*exp(x)+5*x+2;")
%  ylabel('количество итераций')
%  xlabel('Eps')
%  plot(Eps,S_iter)
%  loglog(Eps,S_iter)
% % % % % % % % % 
