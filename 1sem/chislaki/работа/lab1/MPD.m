a1=0;
b1=0.5;
a2=-0.5;
b2=0;
eps=10^(-3);
func1=@(x) x.^4+4*x.^3-12*x.^2+1;
func2=@(x) 5.^x+3*x;
Eps=[];
Roots2=[];
Roots1=[];
Iters1=[];
Iters2=[];
Iters11=[];
Iters22=[];
Approx11=[];
Approx22=[];
for i=1:16
    eps=10^(-1*i);
    Eps=[Eps,eps];
    [koren1,iter1]=MPDq(func1,a1,b1,eps);
    [koren2,iter2]=MPDq(func2,a2,b2,eps);
    Roots1=[Roots1,koren1];
    Roots2=[Roots2,koren2];
    Iters1=[Iters1,iter1];
    Iters2=[Iters2,iter2];
    [koren11,iter11,approx11]=MPDq(func1,a1+i*0.01,b1+i*0.01,10^(-4));
    [koren22,iter22,approx22]=MPDq(func2,a2-i*0.01,b2-i*0.01,10^(-4));
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
legend('полином','трансцендентная')

figure
plot(Approx11,Iters11)
hold on
plot(Approx22,Iters22)
title('Числа итераций от начального приближения')
legend('полином','трансцендентная')

    Approx11=[Approx11,approx11]
    Approx22=[Approx22,approx22]

figure
semilogx(Eps,Iters1)
hold on
semilogx(Eps,Iters2)
title('Числа итераций от заданной точности')
legend('полином','трансцендентная')

figure
hold on
stap=0.00001;
X=[-1:stap:1];
plot(X,func1(X))
plot(koren1,0,'c*')
plot(X,func2(X))
plot(koren2,0,'c*')

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
