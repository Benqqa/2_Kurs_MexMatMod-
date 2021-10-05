function Hess = Hess2(A)
n = length(A);
disp('делаю Хессенберг')
A_ist=A
B=A;
for i=2:n-1
        for j=i+1:n
        c = A(i,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^0.5;
        s = A(j,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^0.5;
        B(i,i)=c^2*A(i,i)+s^2*A(j,j)+2*c*s*A(i,j);
        B(j,j)=s^2*A(i,i)+c^2*A(j,j)-2*c*s*A(i,j);
        B(j,i)=0;
        B(i,j)= B(j,i);
            for m=1:n
                if (m~=i) && (m~=j)
                    B(i,m)=c*A(m,i)+s*A(m,j);
                    B(m,i)=B(i,m);
                    B(j,m)=-s*A(m,i)+c*A(m,j);
                    B(m,j)=B(j,m);
                end
            end    
        end
end
Hess=hess(A)
end