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
    for e=1:15
        eps= 10^(-1*e)
        Eps=[Eps,eps]
        Err=[];
        Iter=[];
        k=1
        D=diag(diag(A));
        E=eye(size(A));
        D1= inv(D);
        C=E-D1*A
        if(norm(C,1)>=1)
        error('1 Ќорма C должна быть <1')
        end
        C_top=zeros(n)
        for i=1:n
            for j= i:n
                   C_top(i,j)=C(i,j)
            end
        end
        C_bott=C-C_top
         Xin=b;
         B=inv(E-C_bott)*C_top
         F=inv(E-C_bott)*D1*b
         u=(1-norm(C,1)/norm(C,1)*eps
        while  k<10000
            Err=[Err,norm(Xin-root)];
            Iter=[Iter,k];
            Xout=B*Xin+F;
            if(norm(Xin-Xout,1)<=(u)break,end 
            k=k+1;
            Xin=Xout;
            
        end
        ChisloIter_eps=[ChisloIter_eps,k];
        Osch=[Osch,norm(Xin-root)/norm(Xin)]
    end
    if q
            ChisloIter=[ChisloIter,k];
    end
end

figure
loglog(DET,ChisloIter,'--*')
xlabel('DET')
ylabel('ChisloIter')
length(Iter)
length(Err)

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