function [lm, ev] = powerrrr(A, xp, L1, L2)
    n = length(A);

    xn = A*xp;
    np = norm(xp, inf);
    xp1 = xp/np;
    xn1 = A*xp1;
    nn = norm(xn1,inf);

    Eps = [];
    Err = [];
    err = abs(np-nn);
    K = [];
    K1 = [];

    k = 0;
    k1 = 0;

    Err = [Err, L2-L1];
    for j = 3:15
        eps = 10^-j;
        while max(Err)>eps
            xp = xn;
            xn = A*xp;
            Err = [];
            for i = 1:n
                L1(i) = L2(i); %без нормировки
                L2(i) = xn(i)/xp(i);
                Err = [Err, L2(i)-L1(i)];
            end
            L1(1)
            k = k+1;
        end

        while err>eps
            xn1 = A*xp1;
            [row, col] = find(abs(xn1) == norm(xn1,inf));
            nn = xn1(row, col);
            xp1 = xn1/nn;
            err = abs(np-nn);
            np = nn; %с нормировкой
            k1 = k1+1;
        end
        Eps = [Eps, eps];
        K = [K, k];
        K1 = [K1, k1];
    end
    lm = nn;
    ev = xp1;
    % emax = [nn, 0, xn1.'];
    figure
    semilogx(Eps, K, '-r*')
    hold on
    grid on
    semilogx(Eps,K1, '-b*')
    title('Зависимость кол-ва итераций от точности')
    legend('без нормировки', 'с нормировкой')
    xlabel('Eps')
    ylabel('Кол-во итераций')
end
