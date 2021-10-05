n=6
A =rand(n,n)
% B =rand(n,n)
% [Q,r]=qr(B);
% d=[1:n];
% D=diag(d);
% A=Q*D*D'*Q'
%%для6 несеметрич
s="не симметричной матрицы"
main(A,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%% симетр
% A =rand(n,n)
% A=(A+A')*1/2
B =rand(n,n)
[Q,r]=qr(B);
d=[0.03 0.06 0.12 0.15 0.21 0.24];
D=diag(d);
A=Q'*D*Q
s="симметричной матрицы"
main(A,s)
%%%%%%%%%%%%%%%%%%% плохо определенные корни
B =rand(n,n)
[Q,r]=qr(B);
d=[1:n-1,n-1+0.000000000000000000001];
D=diag(d);
A=Q'*D*Q
s="мастрицы с плохо отделимыми корнями"
main(A,s)
 axis([-2, 6, -1, 1])

syms x;
ff=det(A-x*eye(n))
GGGG=roots(poly(A))
figure
title('матлаб')
hold on
x6 = -10:0.05:35;
axis([-0, 11, -1, 1])
plot(x6, double(subs(ff, x6)))
plot(GGGG,0,'o')

function m1=main(A,s)
    n = length(A);
    p=LeVe(A)
    syms z y;
    f=p(n+1)
    for i=1:n
        f=f+p(n+1-i)*z^i
    end
    X=NM(p,f,A)
    hold on
    x = -10:0.05:35;

        axis([0, 12, -1, 1])

    figure
    hold on
    grid on
    plot(x, double(subs(f, x)))
    title('применение метода для '+s)
    plot(X,0, '*')
    
end
function X=NM(p,f,A)
    X=[]
    R=[]
    o=1
    u=1
    L_g=[]
    R_g=[]
    [R_g,L_g]=Granizi(f,p,A)
%     R1=roots(p);
%     I=imag(R1);
%     for i=1:length(R1)
%         if I(i)==0
%             R(o)=R1(i) 
%             o=o+1
%         end
%     end
%     for i=1:length(R)
%         R_g(i)=R(i)+0.002
%         L_g(i)=R(i)-0.002
%     end
    for i=1:length(R_g)
        a=L_g(i)
        b=R_g(i)
        fun=inline(f);
        if fun(a)*fun(b)<0 
            eps=10^(-5);
            n_iter=0;
            fun_d=diff(f);
            fun_d2=diff(fun_d);
            
            f_d=inline(fun_d);
            f_d2=inline(fun_d2);
            Y=fun(X);
            Y_d=f_d(X);
            Y_d_2=f_d2(X);
            x0=b;
            m=1
            stap=10^(-7);
            while fun(x0)*f_d2(x0)<=0 || o<2000
                o=o+1;
                x0=x0-stap;
                if o>2000
                    break
                end
            end
            xi=x0;
            xi_1=0;
            k=0;
            disp("корень есть на ["+a+","+b+"]")
            if ((Y_d>0) | (Y_d*-1>0))
                disp("первая производная знакопостяонна на ["+a+","+b+"]")
                if ((Y_d_2>0) | (Y_d_2*-1>0))
                    disp("вторая производная знакопостяонна на ["+a+","+b+"]")
                    while abs(xi-xi_1)>=eps || k==100
                        k=k+1;
                        xi_1=xi;
                        xi=xi-m*(fun(xi)/f_d(xi));
                        if k>1000
                            break
                        end
                    end
                    disp(k)
                end
            end
            X(u)=xi
            u=u+1
        end
    end
end
function p=LeVe(A)
    n = length(A);
    S=[]
    Ai=1
    for i=1:n
        Ai=Ai*A
        S=[S,-trace(Ai)]
    end
    p=ones(1,n+1);
    p(2)=S(1)
    for k=3:n+1
        g=0
        for i=2:k-1
           g=g+p(i)*S(k-i) 
        end
        p(k)=(1/(k-1))*(S(k-1)+g)
    end
end
function [R_g,L_g]=Granizi(f,p,A)
    n=length(p)
    R_g=[]
    L_g=[]
    fun=inline(f);
    MaxG=MAxGran(p)
    p2=p
    if mod(n,2)==0
        for i=2:2:n-1
            p2(i)=p(i)*(-1)
        end
    else
        for i=1:2:n-1
            p2(i)=p(i)*(-1)
        end
    end
    MinG=MAxGran(p2)
%     MinG=0
    ip=MinG;
    k=0;
    %%%%
%     N_A=norm(A)
%     MinG=-1*N_A
%     MaxG=N_A
    %%%%
    step=0.091;
    for i=MinG+step:step:MaxG
        if fun(i)*fun(ip)<0
          k=k+1;
          R_g(k)=i;
          L_g(k)=ip;
        end
        ip=i;
    end
end
function MaxG=MAxGran(p1)
    n=length(p1)
    a0=p1(1)
    m=0
    maxa=0
    for i=2:n-1
        if p1(i)<0
            if  m==0 
                m=i
            end
            if abs(p1(i))>maxa
                maxa=abs(p1(i))
            end
        end
        
    end
    MaxG=1+(maxa/a0)^(1/m)
end






