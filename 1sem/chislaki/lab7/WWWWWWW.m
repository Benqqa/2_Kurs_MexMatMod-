n=5
A=rand(n,n)
d = [n:-1:1];
% d = [(n-1)+0.000000000000000001,n-1:-1:1];
d=[4,3,2,1,i]
A=generate_m(n,d);
A_ist=A;
figure
LU(A_ist,'случайная матрица с заданными с.ч')
LU(Hessenberg(A_ist),'Хессенберг')
legend('матрица','Хессенберг')
function [L,R]=Lr(A)
        n = length(A);
        L = eye(n);
        R = eye(n);
        for m = 1:n
            sm1 = 0;
            sm2 = 0;
            for j = m:n 
                for k = 1:m-1
                    sm1 = sm1+L(m, k)*R(k, j);
                end
                R(m, j) = A(m, j)-sm1;
                sm1 = 0;
            end

            for i = m+1:n
                for k = 1:m-1
                    sm2 = sm2+L(i, k)*R(k, m);
                end
                L(i, m) = (A(i, m)-sm2)/R(m,m);
                sm2 = 0;
            end
        end 
end
function H = Hessenberg(A)
n = length(A);
    for i=2:n-1
            for j=i+1:n
                c = A(i,i-1)/sqrt((A(i,i-1))^2 + (A(j,i-1)^2));
                s = A(j,i-1)/sqrt((A(i,i-1))^2 + (A(j,i-1)^2));
                T = eye(n);
                T(i,i) = c;
                T(j,j) = c;
                T(i,j)= -s;
                T(j,i) = s;
                A = T'*A*T;        
            end
    end
H=A
end
function Lu=LU(A,tit)
Eps = [];
Iter = [];
n = length(A);
A_ist=A
for t = 3:15
        eps = 10^-t;
        k = 0;
        A = A_ist;
        [L0,R0]=Lr(A_ist)
        A_ist=R0*L0
        A=A_ist
        dd1=diag(A_ist);
    while max(tril(A,-1),[],'all')>eps 
        [L,R] =Lr(A)
        A=R*L
        dd2=diag(A);
        dA = abs(abs(dd1(2)) - abs(dd2(2))); 
%         if dA<eps
%             break
%         end
        k = k+1
    end
    Eps = [Eps, eps];
    Iter = [Iter,k];
end
semilogx(Eps,Iter)
hold on
grid on
end
function nA=generate_m(n,d)
    A=rand(n,n);
    [Q,R]=qr(A);
    D=diag(d);
    A=Q*D*Q'
        if det(A)~=0
            m=1;
            for k=0:(n-1)
                M=A(1:k,1:k);
                m=m*det(M);
            end
            if m~=0
                nA=A;        
            else
                generate_m(n,d)
            end
        else
                generate_m(n,d)
        end
end