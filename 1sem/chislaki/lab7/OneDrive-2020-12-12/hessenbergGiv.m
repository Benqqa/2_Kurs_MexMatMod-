n=5
A=rand(n,n)


h=hess(A)
H1=hessenbergGiv1(A)
function Hess = hessenbergGiv1(A)
n = length(A);
% Q = eye(n);

for i = 2:n-1
    for j = i+1:n
%         B=A
        c = A(i,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^0.5;
        s = A(j,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^0.5;
        A(j,i-1)=0;
%         A
%         q = eye(n);
%         q(i,i) = c;
%         q(j,j) = c;
%         q(i,j)= -s;
%         q(j,i) = s;
%         A = q.'*A*q;
%         Q = Q*q;
            for k=1:n
                for l=1:n
                    if k~=j || k~=i || l~=j || l~=i
                        B(k,l)=A(k,l);
                    end
                    if k==i || l==i
                        B(k,l)=A(k,i)*c+A(k,j)*s;
                    end
                    if k==j || l==j
                        B(k,l)=-A(k,i)*s+A(k,j)*c;
                    end
                end
            end
            A(i,i-1)=0;
            for k=1:n
                for l=1:n
                    if k~=j || k~=i || l~=j || l~=i
                       A(k,l)=B(k,l);
                    end
                    if k==i || l==i
                        A(k,l)=B(i,l)*c+B(j,l)*s;
                    end
                    if k==j || l==j
                        A(k,l)=-B(i,l)*s+A(j,l)*c;
                    end
                end
            end
    end
end
Hess = A;
end
