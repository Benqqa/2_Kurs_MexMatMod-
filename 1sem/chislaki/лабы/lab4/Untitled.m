
A = randn(n,n);
A = 0.5*(A+A'); 
A = A + n*eye(n)
A_ist=A
    root=ones(n,1);
    b=A*root
            eps= 10^(-10)
N=size(A,1);
% приведение системы к нормальному виду

C=A'*A;

D=A'*b;

% приведение системы к виду пригодному для итерационного процесса

% Зейделя

for i=1:N

D1(i)=D(i)/C(i,i);

end;

D1=D1'; % транспонирование матрицы

d1=D1;

for i=1:N

for j=1:N

if i==j

C1(i,j)=0;

else

C1(i,j)=-C(i,j)/C(i,i);

end;

end;

end;

% решение системы линейных уравнений методом Зейделя

R1=d1;

while Flag==0

for i=1:N

v=C1(i,1:N); % выделение i-oй строки матрицы

a=dot(v',d1); % вычисление скалярного произведения

d1(i)=a+D1(i);

end;

R2=d1;

s=max(abs(R2-R1));

if s<eps

z1=d1;

z2=s;

return

end;

R1=R2;

end;