n=6
A =rand(n,n);
A=(A+A')*1/2
eps=10^(-5);
A_ist=A

D=diag(diag(A))
B=A-D
% % S1=0
% % S2=0
% % for i=1:n
% %     S1=S1+(D(i,i))^2;
% %     for j=1:n
% %         S2=S2+B(i,j);
% %     end
% % end

% %%выбор оптимального
% r=zeros(n,1)
% for i=1:n
%     for j=1:n
%         r(i)=r(i)+(B(i,j))^2;
%     end
% end
% A
% [rmax,l]=max(r)
% [xmax,k]=max(B(l,:))
% u=2*A(k,l)/(A(k,k)-A(l,l))
% al=sqrt(1/2*(1+1/sqrt(1+u^2)))
% bet=sign(u)*sqrt(1/2*(1-1/sqrt(1+u^2)))
o=0
while o~=100
    o=o+1
    D=diag(diag(A))
    B=A-D
    S1=0
    S2=0
    for i=1:n
        S1=S1+(D(i,i))^2;
        for j=1:n
            S2=S2+B(i,j);
        end
    end
    r=zeros(n,1)
    for i=1:n
        for j=1:n
            r(i)=r(i)+(B(i,j))^2;
        end
    end
    [rmax,idy]=max(r)
    [xmax,idx]=max(B(idy,:))
    while abs(A(i,j))>eps
     for i=1:n-1
        for j=1:n
            u=2*A(i,j)/(A(i,i)-A(j,j))
            c=sqrt(1/2*(1+1/sqrt(1+u^2)))
            s=sign(u)*sqrt(1/2*(1-1/sqrt(1+u^2)))
%             A(i,j)=(A(i,i)-A(j,j))*al*bet/(al^2-bet^2)
            for k=i:n
                t=A(i,k);
                A(i,k)=(c*A(i,k)+s*A(j,k));
                A(j,k)=(-s*t+c*A(j,k))
                eig(A_ist)
            end
        end
     end
    end
end
A_ist
A
