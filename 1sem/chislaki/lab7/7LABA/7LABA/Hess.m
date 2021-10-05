function B = Hess(n, A)
B = zeros(n);
v = zeros(n,1);
P = zeros(n);

for k = 1 : n - 1
   s = -sign(A(k+1,k))*sqrt(sum((A(k+1:n,k)).^2));
   r = sqrt(0.5*(s*s - A(k+1,k)*s));
   v(k) = 0;
   v(k+1) = (A(k+1,k)-s)/(2*r);
        for j = k + 2 : n
        v(j) = A(j,k)/(2*r);
        end
   P = eye(n) - 2*v*v';
   B = P*A*P';
   for m = 1 : n-2
   A = B;
   end
end

