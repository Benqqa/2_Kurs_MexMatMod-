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

Roots2p=[];
Roots1p=[];
Iters1p=[];
Iters2p=[];
Approx11p=[];
Approx22p=[];
Iters11p=[];
Iters22p=[];
for i=1:15
    eps=10^(-1*i);
    Eps=[Eps,eps];
    [koren1,iter1, approx1]=MNq(func11,a1,b1,eps);
    [koren2,iter2, approx2]=MNq(func22,a2,b2,eps);
    Roots1=[Roots1, subs(func11,x,koren1)];
    Roots2=[Roots2, subs(func22,x,koren2)];
    Iters1=[Iters1, iter1];
    Iters2=[Iters2, iter2];
    [koren11,iter11,approx11]=MNq(func11,a1,0.31+i*0.02,10^(-4));
    [koren22,iter22,approx22]=MNq(func22,a2,0.24+i*0.04,10^(-4));
    Iters11=[Iters11,iter11];
    Iters22=[Iters22,iter22];
    Approx11=[Approx11,approx11]
    Approx22=[Approx22,approx22]
    %%%
    [koren1,iter1]=MPDq(func1,a1,b1,eps);
    [koren2,iter2]=MPDq(func2,a2,b2,eps);
    Roots1p=[Roots1p,func1(koren1)];
    Roots2p=[Roots2p,func2(koren2)];
    Iters1p=[Iters1p,iter1];
    Iters2p=[Iters2p,iter2];
    [koren11,iter11,approx11]=MPDq(func1,a1,0.31+i*0.02,10^(-4));
    [koren22,iter22,approx22]=MPDq(func2,a2,0.24+i*0.04,10^(-4));
    Iters11p=[Iters11p,iter11];
    Iters22p=[Iters22p,iter22];
    Approx11p=[Approx11p,approx11]
    Approx22p=[Approx22p,approx22]
end

figure
Eps
Roots1
Roots2
Roots1p
Roots2p
loglog(Eps,Roots1,'*')
hold on
loglog(Eps,Roots2,'*')
loglog(Eps,Roots1p, '*')
loglog(Eps,Roots2p, '*')
title('значение от заданной точности')
legend('полином NM','трансцендентна€ NM','полином MPD','трансцендентна€ MPD')

figure
plot(Approx11,Iters11)
hold on
plot(Approx22,Iters22)
plot(Approx11p,Iters11p)
plot(Approx22p,Iters22p)

title('„исла итераций от начального приближени€')
legend('полином NM','трансцендентна€ NM','полином MPD','трансцендентна€ MPD')

    Approx11=[Approx11,approx11]
    Approx22=[Approx22,approx22]

figure
semilogx(Eps,Iters1)
hold on
semilogx(Eps,Iters2)
semilogx(Eps,Iters1p)
semilogx(Eps,Iters2p)

title('„исла итераций от заданной точности')
legend('полином NM','трансцендентна€ NM','полином MPD','трансцендентна€ MPD')

figure
hold on
stap=0.00001;
X=[-1:stap:1];
plot(X,func1(X))
plot(koren1,0,'c*')
plot(koren2,0,'c*')
plot(X,func2(X))
legend('полином','трансцендентна€')

function [koren, iter, approx ]= MNq(fun,a,b,eps)
    
    syms x;
    iter=0;
    [M2,m1,f,f_d,f_d2,Y_d,Y_d_2]=differ(fun,a,b);
    
    if subs(fun, x,a)*subs(fun, x,b)<0
    disp("корень есть на ["+a+","+b+"]")
    if ((Y_d>0) | (Y_d*-1>0))
    disp("перва€ производна€ знакопост€онна на ["+a+","+b+"]")
    if ((Y_d_2>0) | (Y_d_2*-1>0))
        disp("втора€ производна€ знакопост€онна на ["+a+","+b+"]")
         x0=b;
         stap=10^(-7);
         %первое приближение
         while f(x0)*f_d2(x0)<=0
                x0=x0-stap;
         end
             approx=x0   
         xi=x0;
         xi_1=0;
         while  0.5*M2/m1*abs(xi-xi_1)^2>=eps
            iter=iter+1
            xi_1=xi;
%             xi=xi-(subs(fun, x,xi)/subs(fun_d, x,xi));
            xi=xi-(f(xi)/f_d(xi));

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
function [M2,m1,f,f_d,f_d2,Y_d,Y_d_2]= differ(fun,a,b)
    syms x;
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
end
function [koren, iter, approx ]= MPDq(fun,a,b,eps)
    iter=0;
    if fun(a)*fun(b)<0
        approx=(a+b)/2;
        while (abs(b-a)/2)>eps
            iter=iter+1;
            c=(a+b)/2;
        if fun(a)*fun(c)<0
            b=c;
        else
            a=c;
        end
        end
    else
        disp('Error with a/b')
    end
    koren=(a+b)/2;
end
function y = fucn1(x)
    y=x.^4+4*x.^3-12*x.^2+1;
end
function y = fucn2(x)
    y=5.^x+3*x;
end
function y = dfunc1(x)
    y=4*x*(-6 + 3*x + x.^2);
end
function y = dfunc2(x)
    y=3 + 5.^x *log(5);
end