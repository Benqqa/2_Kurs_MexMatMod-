a1=0.1;
b1=0.5;
a2=-0.5;
b2=0.1;
func1=@(x) x.^4+4*x.^3-12*x.^2+1;
func2=@(x) 5.^x+3*x;
syms x;
func11=x^4+4*x^3-12*x^2+1;
func22=5^x+3*x;
Eps=[];
Roots2=[];
Roots1=[];
Iters1=[];
Iters2=[];
Approx11=[];
Approx22=[];
Iters11=[];
Iters22=[];
for i=1:15
    eps=10^(-1*i);
    Eps=[Eps,eps];
    [koren1,iter1, approx1]=MNq(func11,a1,b1,eps);
    [koren2,iter2, approx2]=MNq(func22,a2,b2,eps);
    Roots1=[Roots1, koren1];
    Roots2=[Roots2, koren2];
    Iters1=[Iters1, iter1];
    Iters2=[Iters2, iter2];
    [koren11,iter11,approx11]=MNq(func11,a1+i*0.01,b1+i*0.01,10^(-4));
    [koren22,iter22,approx22]=MNq(func22,a2-i*0.01,b2-i*0.01,10^(-4));
    Iters11=[Iters11,iter11];
    Iters22=[Iters22,iter22];
    Approx11=[Approx11,approx11]
    Approx22=[Approx22,approx22]
end

figure
semilogx(Eps,Roots1)
hold on
semilogx(Eps,Roots2)
title('значение от заданной точности')
legend('полином','трансцендентна€')

figure
plot(Approx11,Iters11)
hold on
plot(Approx22,Iters22)
title('„исла итераций от начального приближени€')
legend('полином','трансцендентна€')

    Approx11=[Approx11,approx11]
    Approx22=[Approx22,approx22]

figure
semilogx(Eps,Iters1)
hold on
semilogx(Eps,Iters2)
title('„исла итераций от заданной точности')
legend('полином','трансцендентна€')

figure
hold on
stap=0.00001;
X=[-1:stap:1];
plot(X,func1(X))
plot(koren1,0,'c*')
plot(X,func2(X))
plot(koren2,0,'c*')

function [koren, iter, approx ]= MNq(fun,a,b,eps)
    approx=(a+b)/2;
    syms x;
    iter=0;
    fun_d=diff(fun);
    fun_d2=diff(fun_d);
    f=inline(fun);
    f_d=inline(fun_d);
    f_d2=inline(fun_d2);
    
    stap=10^(-7);
    X=[a:stap:b];
    Y=f(X);
    Y_d=f_d(X);
    Y_d_2=f_d2(X);
    M2=max(abs(Y_d_2))
    m1=min(abs(Y_d))
    
    if subs(fun, x,a)*subs(fun, x,b)<0
    disp("корень есть на ["+a+","+b+"]")
    if ((Y_d>0) | (Y_d*-1>0))
    disp("перва€ производна€ знакопост€онна на ["+a+","+b+"]")
    if ((Y_d_2>0) | (Y_d_2*-1>0))
        disp("втора€ производна€ знакопост€онна на ["+a+","+b+"]")
         x0=b;
         %первое приближение
         while f(x0)*f_d2(x0)<=0
                x0=x0-stap;
         end
         x0=(a+b)/2
        
         xi=x0;
         xi_1=0;
         while  abs(xi-xi_1)>eps
            iter=iter+1
            xi_1=xi;
            xi=double(xi-(subs(fun, x,xi)/subs(fun_d, x,xi)));
            xi;
         end
         koren=xi;
    else
        disp("перва€ производна€ Ќ≈ знакопост€онна на ["+a+","+b+"]")
    end
    else
         disp("втора€ производна€ Ќ≈ знакопост€онна на ["+a+","+b+"]")
    end
    else
        disp('Error with a/b')
    end

end
