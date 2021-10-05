clear, clc

Con=[]
Err2=[]
n=10
% for i_i=1:16
Eps=[]
Iter=[]
A = rand(n,n);
[Q,R]=qr(A);
A=R+R'
Con=[Con,cond(A)];
root=2*ones(n,1);
b=A*root


diagA = diag(A);
A = A-diag(diag(A))

for i=1:14
Err=[]
% A = A*A';
% A = A + n*eye(n)
eps = 10^(-1*i);
Eps=[Eps,eps]
k=0;
x=2*zeros(n,1);
x0 = ones(n,1);
err = 10;
while err > (1-norm(diagA))*eps/norm(diagA)
  k=k+1
  x = (b-A*x0)./diagA;
  Err=[Err,norm(x-root)]
  err = norm(x-x0);
  
    x
    root
  x0 = x;
end

Iter=[Iter,(1-norm(diagA))*eps/norm(diagA)]

end
Err2=[Err2,norm(x-root)/norm(x)]
% end
Con
figure
loglog(1:1:k,Err,'g--*')
xlabel('номер итерации')
ylabel('Err')

figure
loglog(-Eps,Iter,'r--o')
xlabel('Eps')
ylabel('аопостериорная оценка')

% figure
% loglog(Con,Err2,'c--*')
% xlabel('Con')
% ylabel('Err2')