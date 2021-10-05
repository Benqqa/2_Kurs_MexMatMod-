 n=10;

DET=[]
ChisloIter=[]


A = randn(n,n);

A = A + n*eye(n)
A_ist=A
NORMC=[]
for d= -10:1
    A=A_ist*10^(-1*d);
    q=false
   
    root=ones(n,1);
    b=A*root
    Eps=[]
    ChisloIter_eps=[]
    Osch=[]
    k=1
    
    D=diag(diag(A));
    E=eye(size(A));
    D1= inv(D);
    C=E-D1*A
    
     if det(A)~=0
        q=true
        DET=[DET,det(A)]
        NORMC=[NORMC,norm(C,1)];
    end
    for e=1:14
        eps= 10^(-1*e)
        Eps=[Eps,eps]
        Err=[];
        Iter=[];
        k=1

        if(norm(C,1)>=1)
        error('1 Ќорма C должна быть <1')
        end
        x_new=zeros(n,1)
        u=(1-norm(C,1))/norm(C,1)*eps;
        while  k<100
            Err=[Err,norm(x_new-root)];
            Iter=[Iter,k];
            x_old=x_new;
            for i=1:n
                s=0;
                for j=1:n
                    if j~=i
                        s=s-A(i,j)/A(i,i)*x_new(j);
                    end
                end
                x_new(i)=s+b(i)*1/A(i,i);
            end
            if(norm(x_old-x_new,1)<=(u))
                break
            end
            k=k+1
            x_old=x_new;
        end
        ChisloIter_eps=[ChisloIter_eps,k];
        Osch=[Osch,norm(x_new-root)/norm(x_new)]
    end
    if q
            ChisloIter=[ChisloIter,k];
    end
end

figure
loglog(DET,ChisloIter,'--*')
xlabel('DET')
ylabel('ChisloIter')

figure
loglog(DET,NORMC,'--*')
xlabel('DET')
ylabel('norm(C)')

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