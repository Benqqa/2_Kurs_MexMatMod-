a=0.4;
b=0.8;
shag=0.00001
X=[a:shag:b]
Y=[fucn(X)];
Y_d=diff(Y)
Y_d_2=diff(Y,2);
M2=max(abs(Y_d_2))
m1=min(abs(Y_d))
eps=0.000001
%stop=0.5*(M2/m1)*abs(x(i)-x(i-1))<eps%условие выхода
x0=b;
if fucn(a)*fucn(b)<0%проверка на знакопостояство/переменность
    %поиск x0
    x0=b; 
    while fucn(x0)*Y_d_2((abs(a-x0)/shag)-1)<0
        x0=x0-0.00001
        disp(x0)
    end
    xi=[]
    xi(1)=x0
    i=2;
    k=(x0-(fucn(x0)/Y_d(abs(a-x0)/shag)))
    while 0.5*(M2/m1)*abs(xi(i-1)-(fucn(xi(i-1)/Y_d(abs(a-xi(i-1))/shag))))>eps
        xi(i)=(xi(i-1)-(fucn(xi(i-1)/Y_d(abs(a-xi(i-1))/shag))))
        i=i+1
    end
end
