n=10

DET=[]
ChisloIter=[]
NC=[]
SC=[]
%генераци€ матрицы с диагональным преобладанием
A = randn(n,n);
A = A + n*eye(n)
A_ist=A
 k=1;
for d= -10:1
    A=A_ist*10^(-1*d);
    k=1;
    root=ones(n,1);
    b=A*root
    Eps=[];
    ChisloIter_eps=[]
    Osch=[]
    D=diag(diag(A));
    E=eye(size(A));
    %проверка условий
    if isempty(find(diag(A)~=0))
        error('¬ диагонале матрицы ј найден нулевой элемент!! ѕересмотри данные')
    end
    D1= inv(D);
    B=E-D1*A
    
    G=D1*b
        if(norm(B,1)>=1)
        error('Ќорма B должна быть <1')
    end
     if det(A)~=0 
        DET=[DET,det(A)]
    end
    
    for e=1:15
    eps= 10^(-1*e)
    Eps=[Eps,eps]
    k=1;
    x_new=b
    %метод €коби
    Err=[];
    Iter=[];
    u=(1-norm(B,1))/norm(B,1)*eps
    while (k<100)
        Err=[Err,norm(x_new-root,1)];
        Iter=[Iter,k];
        xold=x_new;
        for i=1:n
            s=0;
            for j=1:n
                if j~=i
                    s=s+A(i,j)*x_new(j);
                end
            end
            x_new(i)=(1/A(i,i))*(b(i)-s);
        end
        if(norm(x_new-xold,1)<= u)
            break
        end 
        k=k+1;
        x_new=xold;
    end
     ChisloIter_eps=[ChisloIter_eps,k];
     Osch=[Osch,norm(x_new-root,1)/norm(x_new,1)]
    end
    if det(A)~=0
        ChisloIter=[ChisloIter,k];
    end
    x_new
end

figure
loglog(DET,ChisloIter,'--*')
xlabel('DET')
ylabel('ChisloIter')

figure
semilogy(Iter,Err,'--*')
xlabel('Iter')
ylabel('Err')

figure
semilogx(Eps,ChisloIter_eps,'--*')
xlabel('Eps')
ylabel('ChisloIter-eps')

figure
loglog(Eps,Osch,'--*')
xlabel('Eps')
ylabel('относительна€ погрешность')

