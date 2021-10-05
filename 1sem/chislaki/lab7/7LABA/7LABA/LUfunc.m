function [L,U] = LUfunc(A)
k = size(A);
L = zeros(k);
U = eye(k);
if det(A) == 0
    p = 'error'
    return
else
    for i = 1:k
      for j = 1:k
          if i>=j
              sm1 = 0;
                for s = 1:j-1
                   sm1 = sm1 + L(i,s)*U(s,j);
                end   
              L(i,j) = A(i,j) - sm1;
          else
                sm2 = 0;
                for s = 1:i-1
                 sm2 = sm2 + L(i,s)*U(s,j);
                end
             U(i,j) = (1/L(i,i))*(A(i,j) - sm2);
          end
        end
    end
end
end