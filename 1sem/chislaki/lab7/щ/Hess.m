function Hess = Hess(A)
n = length(A);
for i=2:n-1
        for j=i+1:n
        c = A(i,i-1)/sqrt((A(i,i-1))^2 + (A(j,i-1)^2));
        s = A(j,i-1)/sqrt((A(i,i-1))^2 + (A(j,i-1)^2));
        G = eye(n);
        G(i,i) = c;
        G(j,j) = c;
        G(i,j)= -s;
        G(j,i) = s;
        A = G'*A*G;        
        end
end
Hess=A
end