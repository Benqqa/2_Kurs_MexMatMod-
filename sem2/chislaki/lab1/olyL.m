fun= @(x) exp(-x.^2)
a1=-3
b1=3
a2=4
b2=5
n=5
ErrR1=[]
ErrR2=[]
ErrC1=[]
ErrC2=[]
%  figure
for i=1:n
    [errR1,errC1]=interpol(a1,b1,fun,i)
    [errR2,errC2]=interpol(a2,b2,fun,i)
    ErrR1=[ErrR1,errR1]
    ErrR2=[ErrR2,errR2]
    ErrC1=[ErrC1,errC1]
    ErrC2=[ErrC2,errC2]
end
figure
hold on
plot(2:n+1,ErrR1,2:n+1,ErrC1)
legend('равномерная сетка', 'чебышевская сетка')
title('max отклонение интерполяции от числа узлов на [-3,3]')
figure
hold on
plot(2:n+1,ErrR2,2:n+1,ErrC2)
legend('равномерная сетка', 'чебышевская сетка')
title('max отклонение интерполяции от числа узлов на [4,5]')
% for i=1:10
%     interpol(a1,b1,fun,i)
% end

% plot(a:0.0001:b,fun(a:0.0001:b))
function [ErrR,ErrC]=interpol(a,b,fun,n)
gridR=create_gridRavn(a,b,n)
gridC=create_gridCheb(a,b,n)

polyR=polyL(fun,gridR,n)
polyC=polyL(fun,gridC,n)

syms z;
ErrR=max(abs(subs(polyR, z, a:0.01:b)-fun(a:0.01:b)))
ErrC=max(abs(subs(polyC, z, a:0.01:b)-fun(a:0.01:b)))

graph(fun,a,b,polyC,gridC)
end
function g= graph(fun,a,b,poly,grid1)
    syms z;
    figure
    hold on
    grid on
    plot(grid1,fun(grid1),'*')
    plot(a:0.01:b, subs(poly, z, a:0.01:b))
    plot(a:0.0001:b,fun(a:0.0001:b))

end
function [Grid]= create_gridRavn(a,b,n)
    Grid=[];
    h=(b-a)/(n);
    for i=0:n
        Grid=[Grid,a+h*i];
    end
end
function Grid=create_gridCheb(a,b,n)
    Grid=[];
    for i=0:n
        Grid=[Grid,0.5*(a+b)+0.5*(b-a)*cos((2*i+1)*pi/(2*(n+1)))];
    end
end
function [pol] = polyL(fun,grid,n)
    syms z;
    syms m;
    pol=0;
    m=1;
    for i=1:n+1
         m=1;
        for j=1:n+1
            if i~=j
                 m=m*(z-grid(j))/(grid(i)-grid(j));
            end
        end
        pol=pol+fun(grid(i))*m;
    end
end
