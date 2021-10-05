N=100
b=2
T=zeros(N+1, N+1)
for k=1:N
    T(1,k)=100;
    T(k,1)=k*cos(k);
end
for i =2:N-1
    for j=1:N-2
        if i==1
            T(i,j)=100;
        else
            if j==1
                T(i,j)=i*cos(i);
            else
                T(i,j+1)=T(i,j)+b*(T(i+1,j)-2*T(i,j)+T(i-1,j));
                T(N-i,j+1)=T(N-i+1,j)+0.01;
            end
        end
    end
end
[x,y]=meshgrid(1:1:101)
figure
surf(x,y,T)