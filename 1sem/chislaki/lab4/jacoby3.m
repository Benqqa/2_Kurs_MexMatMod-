n=10

DET=[]
ChisloIter=[]
%генераци€ матрицы с диагональным преобладанием
A = randn(n,n);
A = A + n*eye(n)
A_ist=A
for d= -10:1
    A=A_ist*10^(-1*d);
    
    if det(A)~=0 
        DET=[DET,det(A)]
    end
    root=ones(n,1);
    b=A*root
    Eps=[];
    ChisloIter_eps=[]
    Osch=[]
%     for e=1:15
%     eps= 10^(-1*e)
    Eps=[Eps,eps]
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
    k=1;
    Xin=b;
    %метод €коби
    Err=[];
    Iter=[];
    u=(1-norm(B,1))/norm(B,1)*eps
    while (k<1000)
        Err=[Err,norm(Xin-root,1)];
        Iter=[Iter,k];
        Xout=B*Xin+G;
        if(norm(Xin-Xout,1)<= u)
            break
        end 
        k=k+1;
        Xin=Xout;
    end
     ChisloIter_eps=[ChisloIter_eps,k];
     Osch=[Osch,norm(Xin-root,1)/norm(Xin,1)]
    end
    if det(A)~=0
        ChisloIter=[ChisloIter,k];
    end
    
    Xout
% end

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

