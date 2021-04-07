f= @(x)exp(-x.^2)
a1=-1
b1=1

a2=0
b2=1
n=3
Gauss_4(f,a1,b1,n)
function I=Gauss_4(f,a,b,n)
    [X,A]=X_leg(f,n,a,b)
    sum1=0
    for i=1:n
        sim1=sum1+(A(i)*f(X(i)))
    end
    I=sum1
end

function [X,A]=X_leg(f,n,a,b)
X=[]
A=[]
P=legendre(n,a:0.1:b)
X=roots(P)
dP=diff(P)
for i=1:n
    A(i)=2/((1-(X(i))^2)*(dP(X(i)))^2)
end
end