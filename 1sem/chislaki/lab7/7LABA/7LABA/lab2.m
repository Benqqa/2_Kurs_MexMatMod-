A = [1 4 1 3; 0 -1 3 -1; 3 1 0 2; 1 -2 5 1]
b = [1; 2; 3; 4]
disp ('Решение системы A = L*U')
[L, U] = LUfanc (A, b)
disp('Решение системы L*y0 = b')
y0 = forwsub(L,b)
disp('Решение системы U*x0 = y0')
x0 = backsub(U,y0)
     
b1 = disturbance_b(b)
[L, U] = LUfanc (A, b1)
disp('Решение системы L*y1 = b1')
y1 = forwsub(L,b1)
disp('Решение системы U*x1 = y1')
x1 = backsub(U,y1)

A1 = disturbance_A (A, b)
[L, U] = LUfanc (A1, b)
disp('Решение системы L*y2 = b')
y2 = forwsub(L,b)
disp('Решение системы U*x2 = y2')
x2 = backsub(U,y2)

