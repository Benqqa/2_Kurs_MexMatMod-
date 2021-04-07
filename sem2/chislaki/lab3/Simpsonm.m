f= @(x)exp(-x.^2)
a1=-1
b1=1

a2=2
b2=4

Eps=[]
Iter1=[]
Iter2=[]

for i=1:15
eps=10^(-i);
Eps=[Eps,eps];
[iter1,Ans1,Err1]=method_sym(eps,f,a1,b1);
Iter1=[Iter1,iter1];

[iter2,Ans2,Err2]=method_sym(eps,f,a2,b2);
Iter2=[Iter2,iter2];
end

figure
loglog(1:length(Ans1),Ans1,'*-')
title('����������� ������� �� �� ������ ��������')
hold on
loglog(1:length(Ans2),Ans2,'*-')
legend('����. �������','�� ����. �������')
grid on

figure
semilogy(1:length(Err1),Err1,'*-')
title('����������� ������ �� ������ ��������')
hold on
semilogy(1:length(Err2),Err2,'*-')
legend('����. �������','�� ����. �������')
grid on

figure
semilogx(Eps,Iter1,'*-')
title('����������� ����� �������� �� �������� �������')
hold on
semilogx(Eps,Iter2,'*-')
legend('����. �������','�� ����. �������')
grid on

function [iter,Ans,Err]=method_sym(eps,f,a,b)
    Ans=[]
    Err=[]
    I_p=0;
    p=-1
    while true
         p=p+1;
         I=sympson_methods(f,a,b,2^p);
         dI=abs(I-I_p)/15;
         if p>0
             Err=[Err,dI]
             Ans=[Ans,I];
             if dI<eps
                 break
             end
         end
         I_p=I;
    end
    iter=p
end
function ansver=sympson_methods(f,a,b,n)
    N=2*n;
    h=(b-a)/N;
    X=GridRavn(a,b,N)
    sum1=0;
    sum2=0;
    for i =3:2:length(X)-1
        length(X)
        sum1=sum1+f(X(i))
    end
    for i =2:2:length(X)-1
        sum2=sum2+f(X(i))
    end
    ansver=h/3*(f(X(1))+2*sum1+4*sum2+f(X(length(X))));
end
function [Grid]= GridRavn(a,b,n);
    Grid=[];
    h1=(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h1*i];
    end
end
% 1) ��������� ��������
% 2) ���������� �������� 2^n
% 3) n ����� ���������