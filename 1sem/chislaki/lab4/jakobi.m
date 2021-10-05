n=10
A = rand(n,n);
% A = 0.5*(A+A');
A = A*A';
A = A + n*eye(n)
root=ones(n,1);
b=A*root
for r=2:10
a=r*0.1
C=eye(n)-a*A
x=a*b;
eps= 10^(-10)
for i = 1:5
    x_pred=x;
    x=C*x+a*b;
    if norm(C)<=1/2
        if norm(x-x_pred,2)<=eps 
            disp("norm(x-x_pred,2)<=eps")
            break
        end
    else
        if norm(x-x_pred,2)<=(1-norm(C))*eps/norm(C) %аопстериорная оценка
            disp("аопстериорная оценка")
            break
        end
    end 
end
x
end

