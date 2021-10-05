function Hess = qqqqq(A)
n = length(A);
disp('делаю Хессенберг')
B=A;
for i=2:n-1
        for j=i+1:n
            c = A(i,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^0.5;
            s = A(j,i-1)/((A(i,i-1))^2 + (A(j,i-1)^2))^0.5;
            for k=1:n
                for l=1:n
                    if k~=i && k~=j && l~=i && l~=j
                        B(k,l)=A(k,l)
                    end
                    if k==i || l==i 
                        B(k,l)=A(k,i)*c+A(k,j)*s
                    end
                     if k==j || l==j 
                        B(k,l)=-A(k,i)*s+A(k,j)*c
                     end
                end
            end
            for k=1:n
                for l=1:n
                    if k~=i && k~=j && l~=i && l~=j
                        A(k,l)=B(k,l)
                    end
                    if k==i || l==i 
                        A(k,l)=B(i,l)*c+B(j,l)*s
                    end
                     if k==j || l==j 
                        A(k,l)=-B(i,l)*s+B(j,l)*c
                     end
                end
            end
            
        end
end
% Hess=hess(A)
Hess=A
end