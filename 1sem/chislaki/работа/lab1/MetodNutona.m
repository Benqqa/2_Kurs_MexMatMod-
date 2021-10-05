 syms x;
 fun=x^2-x-2;
 eps=10^-4
   a2=-0.5;
    b2=0.1;
    iter=0;
    fun_d=diff(fun);
    fun_d2=diff(fun_d);
    f=inline(fun);
    f_d=inline(fun_d);
    f_d2=inline(fun_d2);
    
    stap=10^(-7);
    X=[a:stap:b];
    Y=f(X);
    Y_d=f_d(X);
    Y_d_2=f_d2(X);
    M2=max(abs(Y_d_2))
    m1=min(abs(Y_d))
    
    if subs(fun, x,a)*subs(fun, x,b)<0
    disp("корень есть на ["+a+","+b+"]")
    if ((Y_d>0) | (Y_d*-1>0))
    disp("перва€ производна€ знакопост€онна на ["+a+","+b+"]")
    if ((Y_d_2>0) | (Y_d_2*-1>0))
        disp("втора€ производна€ знакопост€онна на ["+a+","+b+"]")
         x0=b;
         %первое приближение
         while f(x0)*f_d2(x0)<=0
                x0=x0-stap;
         end
             approx=x0   
         xi=x0;
         xi_1=0;
         while  0.5*M2/m1*abs(xi-xi_1)^2>=eps
            iter=iter+1
            xi_1=xi;
%             xi=double(xi-(subs(fun, x,xi)/subs(fun_d, x,xi)));
            xi=xi-(f(xi)/f_d(xi));

            xi
         end
         koren=xi;
    else
        disp("перва€ производна€ Ќ≈ знакопост€онна на ["+a+","+b+"]")
    end
    else
         disp("втора€ производна€ Ќ≈ знакопост€онна на ["+a+","+b+"]")
    end
    else
        disp('Error with a/b')
    end

