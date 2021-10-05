n=6
A =rand(n,n)
A=(A+A')*1/2

F=norm(A,'fro')
 eps=10^(-5);
 r=0
 while r~=200
    D=diag(diag(A));
    B=A-D;
    [ymax, yind] = max(B);
    [xmax, xind] = max(ymax);
    yind=yind(xind);
     r=r+1;
    for j=1:n-1
        for k=j+1:n

                if A(j,j)==A(k,k)

                    theta = pi/4;
                    c = cos(theta);
                    s = sin(theta);
                else
                  tg = (2*A(xind,yind))/(A(j,j)-A(k,k)) ;
                  c2 = 1/(1+tg^2) ;
                  c=sqrt(c2);
                  s = sqrt(1-c2);
                end
                
%                 A(xind,yind)=c*A(xind,yind)+1/2*s*(A(k,k)-A(j,j));
                T=eye(n)
                T(xind,xind)= c
                T(yind,yind)= c
                T(xind,yind)= s
                T(yind,xind)= -s
              
        end
    end
 
  while r~=200
        
  end
 A