f= @(x)exp(-x.^2);
a1=-1
b1=1

a2=0
b2=2

n=4
eps=10^(-10)
%[iter,Ans,Err]=method_Gauss4(eps,f,a,b,n)
Eps=[]
Iter1=[]
Iter2=[]

for i=1:5
eps=10^(-i);
Eps=[Eps,eps];
[iter1,Ans1,Err1]=method_Gauss4(eps,f,a1,b1,n)
Iter1=[Iter1,iter1];

[iter2,Ans2,Err2]=method_Gauss4(eps,f,a2,b2,n)
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
function [iter,Ans,Err]=method_Gauss4(eps,f,a,b,n)
    
    Ans=[]
    Err=[]
    I_p=0;
    p=-1
    while true
         p=p+1;
         I=Gauss4_method(f,a,b,2^p,n);
         dI=abs(I-I_p);
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
function [I]=Gauss4_method(f,a,b,k,n)
    I=0   
    X=GridRavn(a,b,k)
    sum=0
    for i=1:length(X)-1
        sum=sum+Gauss4(f,X(i),X(i+1),n)
    end
    disp('sum')
    I=sum
end
function [I]=Gauss4(f,a,b,n)
syms x;
%�� ��������� ������ ��� ������� �����������
poly=getPoly(a,b,n)
X=zeros(n,1);
%X=getRoots(poly)
X=solve(poly)
A=getFactors(a,b,n,X)
I= double(getSum(f,A,X,n))
end
function [poly]=getPoly(a,b,n)
    syms x;
    poly=(x-a)^n*(b-x)^n;
    poly=expand(poly);
    poly=diff(poly,x, n);
    poly=expand(poly);
end
function [A]=getFactors(a,b,n,X)
syms t a1 b1;
g=zeros(n,1);
for i=1:n
    f_n= t^(i-1);
    k=subs(int(f_n,a1,b1),a1,a);
    g(i)=subs(k,b1,b);
end
M=X'.^(0);
for i=1:n-1
    M=cat(1,M,X'.^(i));
end
A=M\g;
A=A(1:n,1);
end
function [sum]=getSum(f,A,X,n)
    sum=0;
    for i=1:n
        sum=sum+A(i)*f(X(i));
    end
end
function [Grid]= GridRavn(a,b,n);
    Grid=[];
    h1=(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h1*i];
    end
end