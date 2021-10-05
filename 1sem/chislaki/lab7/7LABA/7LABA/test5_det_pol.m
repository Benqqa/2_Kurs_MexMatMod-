n = 5;
%C = generate_matrix(n)
index = 2;

%тестовый пример
% n = 3;
% A = [5 1 -3;
%     3 0 -2;
%     -4 -1 1]

% B = Hess(n, A)
d = eig(C)
A = (C - d(index)*eye(n));
Z = [];
for i = 1:n
    Z(i) = 1e-300;
end
Z = Z';
for i = 1:n
   if A(i, i) == 0
       A(i, i) = 1e-300;
   end
end
A = [A Z];
X = res(A);
X = X ./ norm(X)
