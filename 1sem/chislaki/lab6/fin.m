n=5

d1=[1:n];
d=[1,2,2.000000001,3,1.00000000000001];
dl(d,n,'хорошо')
dl(d1,n,'плохо')
function dl(d,n, tit)
Eps=[]
A =rand(n,n);


A=(A+A')/2
% A=[1.43322109803814,0.00209397012585013,0.133544377714551,0.137678873190195,0.456887973796413;0.00209397012585013,1.99032275382563,0.0522461877523039,-0.0824103374109570,0.00762132030068798;0.133544377714551,0.0522461877523039,1.04385495821882,0.0381120719645597,0.141121953909356;0.137678873190195,-0.0824103374109570,0.0381120719645597,1.05072383183220,0.144746187693895;0.456887973796413,0.00762132030068798,0.141121953909356,0.144746187693895,1.48187735808521]
[Q,R]=qr(A);
D=diag(d);
A=Q*D*Q'
o=0;
A_ist=A
Iter=[]
eps2=10^(-16)
r=zeros(1,n);
    D=diag(diag(A_ist));
    K1=A_ist-D;
    K2=K1
    for i=1:n
        for j=1:n
           r(i)=r(i)+(K1(i,j))^2;
        end
    end
    r_ist=r
[rmax,j]=max(r);
[xmax,i]=max(abs(K2(j,1:n)));
j_ist=j
i_ist=i
for ep=3:16
    eps2=10^(-1*ep);
    o=0;
    rmax=[]
    xmax=[]
    A=A_ist
    B=A_ist
    r=r_ist
    j=j_ist
    i=i_ist
    while o~=2000 || max(r)<eps2 || S2<eps2
        o=o+1
      
        B=A;
        S1=norm(D,'fro');
        S2=norm(K1,'fro');
        
        if max(r)<eps2
%             Iter=[Iter,o];
            break

        end
        if S2<eps2
%             Iter=[Iter,o];
            break

        end
        
        if A(j,j)==A(i,i) 
            theta = pi/4;
            c = cos(theta);
            s = sin(theta);
        else
            tg = (2*A(i,j))/(A(i,i)-A(j,j)) ;%tg(pi/4 *2)=inf
            c2=1/sqrt(1+tg^2)
            c=sqrt((1+c2)/2)
            s=sign(tg)*sqrt((1-c2)/2)
        end
        B(i,i)=c^2*A(i,i)+s^2*A(j,j)+2*c*s*A(i,j);
        B(j,j)=s^2*A(i,i)+c^2*A(j,j)-2*c*s*A(i,j);
        B(j,i)=0;
        B(i,j)= B(j,i);
            for m=1:n
                if (m~=i) && (m~=j)
                    B(i,m)=c*A(m,i)+s*A(m,j);
                    B(m,i)=B(i,m);
                    B(j,m)=-s*A(m,i)+c*A(m,j);
                    B(m,j)=B(j,m);
                end
            end
          D1=diag(diag(A));
          D2=diag(diag(B));
          K1=A-D1;
          K2=B-D2;
          C= (K1==K2)
          r
          for i1=1:n
              if (i1==i) || (i1==j)
                  r(i1)=0;
                  for k=1:n
                    r(i1)=r(i1)+(K2(i1,k))^2;
                  end
              end
          end
          [rmax,j]=max(r);
           [xmax,i]=max(abs(K1(j,1:n)));
        A=B;
%         o=o+1
    end
    Iter=[Iter,o];
    o=0;
    Eps=[Eps, eps2];
end
Eps
% Iter=sort(Iter)

sort(diag(A))
sort(eig(A_ist))
figure
semilogx(Eps,Iter,'c*-')
title(tit)
xlabel('Eps')
ylabel('Iter')
end
