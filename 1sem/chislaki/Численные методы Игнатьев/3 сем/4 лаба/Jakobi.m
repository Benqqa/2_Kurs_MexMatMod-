n=10;
Default_eig = linspace(1,n,n);
B = [0; 0; 0;0; 0; 0;0; 0; 0;0];

for i = 1:1:n
    vec(i) = sin(i);
end

otrazhenia = eye(n);
for j = 1:1:n
    for k = 1:1:n
        otrazhenia(j,k) = (otrazhenia(j,k) - 2 * vec(j) * vec(k)/ norm(vec,2).^2 );
    end
end
figure();
hold all;
grid on;
for i = 0:1:14
    
    MAT = eye(n);
    for k = 1:1:n
        MAT(k,k) = k;
    end
    A = otrazhenia * MAT * otrazhenia;
    
    matrix = fopen('input_matlab.txt', 'w');
    fprintf(matrix,'%d %d\r\n',n, i);
    for j=1:1:n
        for k=1:1:n
            fprintf(matrix,'%.15f ',A(j,k));
        end
        fprintf(matrix,'\r\n');
    end
    
    fclose(matrix);
    
    system('C:\Users\danii\Desktop\ignatiev\C\numerical methods\Lab4_Jakobi\Debug\Lab4_Jakobi');
    
    input_C = fopen('root.txt', 'r');
    eigen_new = fscanf(input_C, '%f %f %f %f %f %f %f %f %f %f');
    fclose(input_C);
        eigen_new=eigen_new';
        
        input_eig = fopen('eigen.txt', 'r');
    
    X=dlmread('eigen.txt');
    fclose(input_eig);
   % eigen_new1=sort(eigen_new);
    
    for k=1:1:round(n/2)
        if (round(eigen_new(k)) ~= k)
          for j =1:1:n
              temp = -X(j,k);
              X(j,k) = -X(j,n-k);
              X(j,n-k) = temp;
          end
          temp = eigen_new(k);
          eigen_new(k) = eigen_new(n-k);
          eigen_new(n-k) = temp;
        end
    end
    
    eig_mat = eye(n);
    for k=1:1:n
        eig_mat(k,k) = eigen_new(k);
    end
 
    disp(X - otrazhenia)
    disp(norm(X - otrazhenia,inf))
    disp(norm(eigen_new - Default_eig, inf))
    
    First_norm = norm(eigen_new - Default_eig, inf);
    Second_norm = norm(A*X - X*eig_mat, inf);
    Third_norm = norm(otrazhenia - X , inf);
    plot(-i,log10(First_norm),'ro')
    plot(-i,log10(Second_norm),'b*')
    plot(-i,log10(Third_norm),'m*')
end
set(gca,'FontSize',16)
title('График зависимости погрешности от точности.');
xlabel('Точность == 10^x');
ylabel('Погрешность == 10^y');
legend('погрешность||X-X*||','погрешность||AX-?X||','погрешность |||?-?*|','location','northwest');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           2      %%55
figure();
hold all;
grid on;
for i = 1:1:15
    
    MAT = eye(n);
    for k = 1:1:n
        MAT(k,k) = k;
    end
    A = otrazhenia * MAT * otrazhenia;
    
    matrix = fopen('input_matlab.txt', 'w');
    fprintf(matrix,'%d %d\r\n',n, i);
    for j=1:1:n
        for k=1:1:n
            fprintf(matrix,'%.15f ',A(j,k));
        end
        fprintf(matrix,'\r\n');
    end
    
    fclose(matrix);
    
    system('C:\Users\danii\Desktop\ignatiev\C\numerical methods\Lab4_Jakobi\Debug\Lab4_Jakobi');
    
    input_C = fopen('iterations.txt', 'r');
    iter = fscanf(input_C, '%d');
    fclose(input_C);
    
    plot(i,(iter),'r*')
end
set(gca,'FontSize',16)
title('График зависимости итераций от точности.');
xlabel('Точность == 10^-x');
ylabel('Итерации == y');
legend('кол-во итераций','location','northwest');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           3     %%55

figure();
hold all;
grid on;
for i = 0:1:6
    
    MAT = eye(n);
    for k = 1:1:n
        MAT(k,k) = k;
    end
    
    MAT(n-1,n-1) = MAT(n-1,n-1) - 10^-i;
    A = otrazhenia * MAT * otrazhenia;
    
    matrix = fopen('input_matlab.txt', 'w');
    fprintf(matrix,'%d %d\r\n',n, 7);
    for j=1:1:n
        for k=1:1:n
            fprintf(matrix,'%.15f ',A(j,k));
        end
        fprintf(matrix,'\r\n');
    end
    
    fclose(matrix);
    
    system('C:\Users\danii\Desktop\ignatiev\C\numerical methods\Lab4_Jakobi\Debug\Lab4_Jakobi');
    
    input_C = fopen('root.txt', 'r');
    eigen_new = fscanf(input_C, '%f %f %f %f %f %f %f %f %f %f');
    fclose(input_C);
        eigen_new=eigen_new';
        
        input_eig = fopen('eigen.txt', 'r');
    
    X=dlmread('eigen.txt');
    fclose(input_eig);
   % eigen_new1=sort(eigen_new);
    
    for k=1:1:round(n/2)
        if (round(eigen_new(k)) ~= k)
          for j =1:1:n
              temp = -X(j,k);
              X(j,k) = -X(j,n-k);
              X(j,n-k) = temp;
          end
          temp = eigen_new(k);
          eigen_new(k) = eigen_new(n-k);
          eigen_new(n-k) = temp;
        end
    end
    
    eig_mat = eye(n);
    for k=1:1:n
        eig_mat(k,k) = eigen_new(k);
    end
 
    disp(X - otrazhenia)
    disp(norm(X - otrazhenia,inf))
    disp(norm(eigen_new - Default_eig, inf))
    
    First_norm = norm(eigen_new - Default_eig, inf);
    Second_norm = norm(A*X - X*eig_mat, inf);
    Third_norm = norm(otrazhenia - X , inf);
    plot(-i,log10(First_norm),'ro')
    plot(-i,log10(Second_norm),'b*')
    plot(-i,log10(Third_norm),'m*')
end
set(gca,'FontSize',16)
title('График зависимости погрешности от приближения.');
xlabel('приближение == 10^x');
ylabel('Погрешность == 10^y');
legend('погрешность||X-X*||','погрешность||AX-?X||','погрешность |||?-?*|','location','northwest');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           4      %%55
figure();
hold all;
grid on;
for i = 1:1:6
    
    MAT = eye(n);
    for k = 1:1:n
        MAT(k,k) = k;
    end
    MAT(n-1,n-1) = MAT(n-1,n-1) - 10^-(i);
    A = otrazhenia * MAT * otrazhenia;
    
    matrix = fopen('input_matlab.txt', 'w');
    fprintf(matrix,'%d %d\r\n',n, 6);
    for j=1:1:n
        for k=1:1:n
            fprintf(matrix,'%.15f ',A(j,k));
        end
        fprintf(matrix,'\r\n');
    end
    
    fclose(matrix);
    
    system('C:\Users\danii\Desktop\ignatiev\C\numerical methods\Lab4_Jakobi\Debug\Lab4_Jakobi');
    
    input_C = fopen('iterations.txt', 'r');
    iter = fscanf(input_C, '%d');
    fclose(input_C);
    
    plot(i,(iter),'r*')
end
set(gca,'FontSize',16)
title('График зависимости итераций от приближения.');
xlabel('Приближение == 10^-x');
ylabel('Итерации == y');
legend('кол-во итераций','location','northwest');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         5

figure();
hold all;
grid on;
for i = 0:1:6
    
    MAT = eye(n);
    for k = 1:1:n
        MAT(k,k) = 1+k* 10^(-i);
    end

    A = otrazhenia * MAT * otrazhenia;
    
    matrix = fopen('input_matlab.txt', 'w');
    fprintf(matrix,'%d %d\r\n',n, 10);
    for j=1:1:n
        for k=1:1:n
            fprintf(matrix,'%.15f ',A(j,k));
        end
        fprintf(matrix,'\r\n');
    end
    
    fclose(matrix);
    
    system('C:\Users\danii\Desktop\ignatiev\C\numerical methods\Lab4_Jakobi\Debug\Lab4_Jakobi');
    
    input_C = fopen('iterations.txt', 'r');
    iter = fscanf(input_C, '%d');
    fclose(input_C);
    
    plot(-i,(iter),'r*')
end
set(gca,'FontSize',16)
title('График зависимости итераций от отделения с.ч. .');
xlabel('Разница между с.ч.  == 10^x');
ylabel('Итерации == y');
legend('кол-во итераций','location','northwest');
disp('end.')