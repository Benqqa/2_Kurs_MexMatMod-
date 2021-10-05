 a=-2.73;
b= -1.59;
eps=10^(-3);
Err=[];
N_iter=[];
n_iter=0;
root=-1.7321;
stap=0.00001;
X=[a:stap:b];
Eps=[];
S_iter=[];
Pogr=[];
if fucn(a)*fucn(b)<0
% while eps>10^(-15)
% eps=eps/10
% Eps=[Eps,eps];
    n_iter=0;
while (abs(b-a)/2)>eps
    Err=[Err,abs((a+b)/2-root)]
        n_iter=n_iter+1;
        N_iter=[N_iter,n_iter];
        disp("приближение"+n_iter+" = "+(a+b)/2);
        c=(a+b)/2;
        if fucn(a)*fucn(c)<0
            b=c;
        else
            a=c;
        end
end
%     S_iter=[S_iter,n_iter]
end
% end
koren=(a+b)/2
% % % % % % % % % 
% figure
%  hold on
%  grid on
%  title(" y=x.^4 - x.^3 - 2*x.^2 + 3*x-3;")
%%%  title(" y=3*exp(x)+5*x+2;")
%  xlabel('Ox')
%  ylabel('Oy')
% plot(X,fucn(X))
% plot(koren,0,'c*')
% % % % % % % % % 

%  hold on
%  grid on
 title(" y=x.^4 - x.^3 - 2*x.^2 + 3*x-3;")
%  title(" y=3*exp(x)+5*x+2;")
 xlabel('номер итераций')
 ylabel('ошибка')
%  plot(N_iter,Err)
 semilogy(N_iter,Err)
% % % % % % % % % 
% figure
%  hold on
%  grid on
%  title(" y=x.^4 - x.^3 - 2*x.^2 + 3*x-3;")
%  ylabel('количество итераций')
%  xlabel('Eps')
%  plot(Eps,S_iter)
%  loglog(Eps,S_iter)
% % % % % % % % % 
