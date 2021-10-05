n=10
Con=[]
Err=[]
err=[]
m=0
d=0
for i_i=1:16
    %%% ��������� �������
    while true
        A=randn(n,n);
        [Q,R]=qr(A);
        for i=1:n
            D(i,i)=10^(16)*i^(-i_i);
        end
        A=Q*D*Q';
        a=A(1,1);
        b=A(2,1);
        c=a/sqrt(a^2+b^2);
        s=b/sqrt(a^2+b^2);
        %%% �������� ������� ������������
        if (-s*a+c*b==0) && (c^2+s^2==1)
            break
        else
            disp("�������������� ������� �� �������� �� ������� ������������, ���������� �����")
        end
    end
    A_ist=A
    Con=[Con,cond(A)];
    root=ones(n,1);
    %%%
    b=A*root% ������ b
    b_ist=b;
    A=[A,b];
    %%% ����� ��������
    for i=1:n-1
        for j=i+1:n
            a=A(i,i);
            b=A(j,i);
            c=a/sqrt(a^2+b^2);
            s=b/sqrt(a^2+b^2);
            for k=i:n+1
                t=A(i,k);
                A(i,k)=(c*A(i,k)+s*A(j,k));
                A(j,k)=(-s*t+c*A(j,k));
            end
        end
    end
    b=A(:,n+1)
    A=A(:,1:n)
    %�������� ��� ������ �����
    x=zeros(n,1);
    x(n)=b(n)/A(n,n);
    for k=n-1:-1:1
        x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
    end
    Err=[Err,norm(x-root)/norm(x)];
    x
    A=A_ist
    b=b_ist
    for k=0:(n-1)
            M=A(1:k,1:k);
            d=d*det(M)
        end
        if d~=0
            break
        end
    %������ ������������=0
            if det(A)~=0
            %���������� � ����������� �������    
            A=[A,b]


            %������ ��� ������ ������
                for i=1:(n)%������� �����
                    for j=(i+1):(n) %������� ��������

            %�������� �� ������� �������
                        if A(i,i)==0
                            disp("���� �� ���������")
                            max=0; %���� �������
                            index_max=0; %������ ����
                            for r=i+1:n %������� �����
                                if abs(A(r,i))>max %���� ��-� ������ ���������
                                    max=abs(A(r,i)); %�� �� ���������� ����������
                                    index_max=r; %�������� ��� ������
                                end
                            end
                            c=A(i,:);%����� i-��� ������ � �� ����� 
                            A
                            A(i,:)=A(r,:);%������ �� � r-��� �������,�.� ������������
                            A(r,:)=c
                        end
            %�������� ���������

                        m=A(j,i)/A(i,i);
                        for k=1:(n+1)%� �������� �, ��� �������� j-��� ������
                            A(j,k)=A(j,k)-m*A(i,k); %����� ���� ������� � �������� �������
                        end
                    end
                end

            %������ ������������!=0
            else
                disp("det!=0") 
            end
            %��������� ����������� �������
            b=A(:,n+1)
            A=A(:,1:n)

            %�������� ��� ������ ������
            x=zeros(n,1);
            x(n)=b(n)/A(n,n);%������ ����� ������ ���������
            for k=n-1:-1:1 %���� ������ ��� ����������
                x(k)=(b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k); %������� � ������� ���������, ������� �� ���������� � �� ���������
            end
            
   err=[err,(norm(x-root))/norm(x)]
end

figure
loglog(Con,Err,'r--o',Con,err,'g--o')
legend('Rotation', "Gauss")
xlabel('Con')
ylabel('Err')

