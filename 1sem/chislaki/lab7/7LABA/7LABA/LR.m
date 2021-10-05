function [A0, iter] = LR(n, eps, A, L, U)
A0 = A;
L0 = L;
U0 = U;
% eps = 1e-18;
dA = eps + 1;
iter = 0;

while(dA > eps)
A0 = L0*U0;
A = U0*L0;
d1 = diag(A0);
d = diag(A);                          
A0 = A;
[L0, U0] = LUfanc (A);
iter = iter + 1;
dA = abs(abs(d1(2)) - abs(d(2))); 
end

end
