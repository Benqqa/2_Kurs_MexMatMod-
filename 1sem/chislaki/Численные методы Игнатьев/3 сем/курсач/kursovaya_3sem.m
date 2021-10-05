X = [0,1,2,3,4,5,6,7,8,9];
X=X';
n=10;
k = 10;
%disp(X);

for i = 1:1:n
    vec(i) = rand(1);
end

otrazhenia = eye(n);
for j = 1:1:n
    for k = 1:1:n
        otrazhenia(j,k) = otrazhenia(j,k) - 2 * vec(j) * vec(k)/ norm(vec,2).^2;
    end
end
MAT=eye(10);
average=[0 0 0 0 0];
iterations_Zeid = [0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
figure();
hold all
grid on
for i = 1:1:14
    sum =0;
    for test = 1:1:11
        for k = 1:1:n
            MAT(k,k) = (rand(1)*10)*(k - 1) / (9) + 1;
        end
        A = otrazhenia * MAT * otrazhenia;

        matrix = fopen('input_matlab.txt', 'w');
        fprintf(matrix,'%d %d\r\n',n, i);
        for j=1:1:10
            for k=1:1:10
                fprintf(matrix,'%.15f ',A(j,k)); %%%%%%%MATRIX
            end
            fprintf(matrix,'\r\n');
        end

        b = MAT*(X);
        for k=1:1:10
            fprintf(matrix,'%.15f ',b(k));     %%%%%%%%B
        end
        fprintf(matrix,'\r\n');
        for k=1:1:10
            fprintf(matrix,'%.15f ',0);        %%%%%%%%root
        end

        fclose(matrix);

        system('C:\Users\danii\Desktop\ignatiev\C\numerical methods\Kursach_CHM\Debug\Kursach_CHM.exe');

        iterations = fopen('iterations.txt', 'r');
        iter_plus = fscanf(iterations, '%d',1);
        iter_mul = fscanf(iterations, '%d',1);
        temp_iter = fscanf(iterations, '%d',1); 
        sum =sum + temp_iter;
        
        fclose(iterations);

        input_C = fopen('output_C.txt', 'r');

        X_LU = fscanf(input_C, '%f',[1 10]);  %%%%%%%%%%LU ROOT

       % X_new = fscanf(input_C, '%f',[1 10]);  %%%%%%%%%1/2
        %  X_new = X_new';

    %     iterat = [iter_mul/2 iter_mul 2*iter_mul];
    %     Norm = [0 0 0];
    %     
    %     
    %     
    %     Norm(1) = norm(X_new - X_LU, inf);
    %     average(1) = average(1)+ log10(Norm(1));
    %     X_new = fscanf(input_C, '%f',[1 10]);  %%%%%%%%%1
    %     Norm(2) = norm(X_new - X_LU, inf);
    %     average(2) = average(2)+  log10(Norm(2));
    %     X_new = fscanf(input_C, '%f',[1 10]);  %%%%%%%%%2
    %     Norm(3) = norm(X_new - X_LU, inf);
    %     average(3) = average(3)+  log10(Norm(3));

        %  hPlot = plot(log10(Norm),iterat,'k');
        % set( hPlot, 'LineWidth', 2 );

        fclose(input_C);
        disp(i)
        disp(temp_iter)
    end
    
    iterations_Zeid(i) = sum/test;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%operations LU
hPlot = plot([15 0],log10([iter_plus iter_plus]),'r-.');
set( hPlot, 'LineWidth', 3 );

hPlot =  plot([15 0],log10([iter_mul/2 iter_mul/2]),'b-.');
set( hPlot, 'LineWidth', 3 );
hPlot =  plot([15 0],log10([iter_mul iter_mul]),'b-.');
set( hPlot, 'LineWidth', 3 );
hPlot =  plot([15 0],log10([iter_mul*2 iter_mul*2]),'b-.');
set( hPlot, 'LineWidth', 3 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hPlot = plot(-average/i,iterat,'g');
hPlot = plot([1 2 3 4 5 6 7 8 9 10 11 12 13 14], log10(iterations_Zeid*100),'g');
set( hPlot, 'LineWidth', 3 );

set(gca,'FontSize',15)
title('График зависимости операций прямого и итерационного методов от точности.');
ylabel('операции == 10^y');
xlabel('точность == 10^-^х');
legend('Операции суммирования (440)','Операции умножения (324)','Операции умножения (648)','Операции умножения (1296)','Кол-во операций для заданной точности м. Зейделя');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
hold all;
C = logspace(1,10,10);
iter = [0 0 0 0 0 0];

cond = 6;
for i = 1:1:cond
    MAT = eye(10);
    for k = 1:1:10
        MAT(k,k) = (k - 1) / (10 - 1) * (C(i) - 1) + 1;
    end
    A = otrazhenia * MAT * otrazhenia;

   matrix = fopen('input_matlab.txt', 'w');
        fprintf(matrix,'%d %d\r\n',n, 15);
        for j=1:1:n
            for k=1:1:n
                fprintf(matrix,'%.15f ',A(j,k));
            end
            fprintf(matrix,'\r\n');
        end

        b = A*X;
        for k=1:1:n
            fprintf(matrix,'%.15f ',b(k));
        end
        fprintf(matrix,'\r\n');
        for k=1:1:n
            fprintf(matrix,'%.15f ',0);
        end
        
        fclose(matrix);
        
        system('C:\Users\danii\Desktop\ignatiev\C\numerical methods\LAB_3_Zeidel\Debug\LAB_3_Zeidel.exe');
    
    input_C = fopen('iterations.txt', 'r');
    
    iter(i) = fscanf(input_C, '%d');
    fclose(input_C);
    disp(num2str(i))
end

hPlot = plot([1 2 3 4 5 6],log10(iter/2),'r');
set( hPlot, 'LineWidth', 2 );

hPlot = plot([1 cond],log10([iter_plus iter_plus]),'g-.');
set( hPlot, 'LineWidth', 2 );

hPlot = plot([1 cond],log10([iter_mul iter_mul]),'b-.');
set( hPlot, 'LineWidth', 2 );

set(gca,'FontSize',10)
title('График зависимости операций м. Зейделя от cond.');
xlabel('cond == 10^x');
ylabel('итерации == 10^y');
legend('Операции умножения м. Зейделя','Операции суммирования (440)','Операции умножения (628)');

disp('END.')
