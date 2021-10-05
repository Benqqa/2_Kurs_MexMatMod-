n=10
A=rand(n,n)
% d = [n:-1:1];
d = [(n-1)+0.000000000000000001,n-1:-1:1];
A=genA(n,d);
A0=A;
figure
LU(A0,'обычная матица')
LU(Hess(A0),'приведеная к Хессенбергу методом вращений')
legend('обычная матица','приведеная к Хессенбергу методом вращений')
function [L,R]=Lr(A)
        n = length(A);
        L = eye(n);
        R = eye(n);
        for m = 1:n
            S1 = 0;
            S2 = 0;
            for j = m:n 
                for k = 1:m-1
                    S1 = S1+L(m, k)*R(k, j);
                end
                R(m, j) = A(m, j)-S1;
                S1 = 0;
            end

            for i = m+1:n
                for k = 1:m-1
                    S2 = S2+L(i, k)*R(k, m);
                end
                L(i, m) = (A(i, m)-S2)/R(m,m);
                S2 = 0;
            end
        end 
end
function Hess1 = Hess(A)
n = length(A);
    for i=2:n-1
            for j=i+1:n
            c = A(i,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^(1/2);
            s = A(j,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^(1/2);
            G = eye(n);
            G(i,i) = c;
            G(j,j) = c;
            G(i,j)= -s;
            G(j,i) = s;
            A = G'*A*G;        
            end
    end
Hess1=A
end
function LUL=LU(A,tit)
EPS = [];
N_Iter = [];
n = length(A);
A0=A
for t = 3:15
        eps = 10^-t;
        k = 0;
        A = A0;
        [L0,R0]=Lr(A0)
        A=A0
        v=L0*R0-A0;
        A0 = R0*L0-v
         A=A0
        dd1=diag(A0);
    while max(tril(A,-1),[],'all')>eps 
        [L,R] =Lr(A)
        v=L*R-A;
        A = R*L-v
        dd2=diag(A);
        dA = abs(abs(dd1(2)) - abs(dd2(2))); 
        if dA<eps
            break
        end
        k = k+1;
    end
    EPS = [EPS, eps];
    N_Iter = [N_Iter,k];
end
semilogx(EPS,N_Iter)
hold on
grid on
end
function new_A=genA(n,d)
    A=rand(n,n);
    [Q,R]=qr(A);
    D=diag(d);
    A=Q*D*Q'
    while true
        if det(A)~=0
            m=1;
            for k=0:(n-1)
                M=A(1:k,1:k);
                m=m*det(M);
            end
            if m~=0
                new_A=A;
                break
            else
                genA(n,d)
                
            end
        else
                genA(n,d)
                
        end
        
    end
    
end