 n=10;

DET=[]
ChisloIter=[]

A = randn(n,n);
A = 0.5*(A+A'); 
A = A + n*eye(n)
A_ist=A
for d= -10:1
    A=A_ist*10^(-1*d);
    q=false
    if det(A)~=0
        q=true
        DET=[DET,det(A)]
    end
    root=ones(n,1);
    b=A*root
    Eps=[]
    ChisloIter_eps=[]
    Osch=[]
     k=1
    for e=1:14
        eps= 10^(-1*e)
        Eps=[Eps,eps]
        Err=[];
        Iter=[];
        k=1
        D=diag(diag(A));
        E=eye(size(A));
        D1= inv(D);
        C=E-D1*A
        C_top=zeros(n)
        for i=1:n
            for j= i:n
                   C_top(i,j)=C(i,j)
            end
        end
        C_bott=C-C_top
        B=(E-C_bott)^(-1)*C_top
        if(norm(C,1)>=1)
        error('1 Ќорма C должна быть <1')
        end
        x_new=zeros(n,1)
        u=(1-norm(B,1))/norm(B,1)*eps;
        while  k<100
            Err=[Err,norm(x_new-root)];
            Iter=[Iter,k];
            x_old=x_new;
            for i=1:n
                s=0;
                for j=i+1:n
                        s=s-A(i,j)/A(i,i)*x_old(j);
                end
                for j=1:i-1
                        s=s-A(i,j)/A(i,i)*x_new(j);
                end
                x_new(i)=s+(1/A(i,i))*b(i);
            end
            if(norm(x_old-x_new,1)<=(u))
                break
            end
            k=k+1
            x_old=x_new;
        end
        ChisloIter_eps=[ChisloIter_eps,k];
        Osch=[Osch,norm(x_new-root,1)/norm(x_new,1)]
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