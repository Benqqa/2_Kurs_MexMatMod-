n=5
A =rand(n,n);
Eps=[]

A=(A+A')/2
A=[1.43322109803814,0.00209397012585013,0.133544377714551,0.137678873190195,0.456887973796413;0.00209397012585013,1.99032275382563,0.0522461877523039,-0.0824103374109570,0.00762132030068798;0.133544377714551,0.0522461877523039,1.04385495821882,0.0381120719645597,0.141121953909356;0.137678873190195,-0.0824103374109570,0.0381120719645597,1.05072383183220,0.144746187693895;0.456887973796413,0.00762132030068798,0.141121953909356,0.144746187693895,1.48187735808521]
% d=[1,2,3,4,5];
d=[1,2,2.000000001,3,1.00000000000001];
[Q,R]=qr(A);
D=diag(d);
A=Q.'*D*Q
A_ist=A

A_ist=A
Iter=[]
for ep=3:14

eps=10^(-ep);
Eps=[eps,Eps]
o=0
% 
% D=diag(diag(A));
%     K=A-D
%     for i=1:n
%         for j=1:n
%            r(i)=r(i)+(K(i,j))^2;
%         end
%     end
while o~=2000
    o=o+1
     B=A;
    D=diag(diag(B));
    K=B-D
    S1=norm(D,'fro');
    S2=norm(K,'fro')
    r=zeros(1,n);
    for i=1:n
        for j=1:n
           r(i)=r(i)+(K(i,j))^2;
        end
    end
    if max(r)<eps
        break
    end
    if S2<eps
        break
    end
    [rmax,j]=max(r);
    [xmax,i]=max(abs(K(j,1:n)));
    if A(j,j)==A(i,i)
        theta = pi/4;
        c = cos(theta);
        s = sin(theta);
        
    else
        tg = (2*A(i,j))/(A(i,i)-A(j,j)) ;
        c2=1/sqrt(1+tg^2)
        c=sqrt((1+c2)/2)
        s=sign(tg)*sqrt((1-c2)/2)
    end

    G=eye(n)
    G(i,j)=-s;
    G(i,i)=c;
    G(j,j)=c;
    G(j,i)=s;
    A=G'*A*G;
end
Iter=[Iter,o]
end
Eps
% A_ist
% A
sort(diag(A))
sort(eig(A_ist))
% Err=sort(diag(A))-sort(eig(A_ist))
figure
semilogx(Eps,Iter,'c*-')
xlabel('Eps')
ylabel('Iter')
