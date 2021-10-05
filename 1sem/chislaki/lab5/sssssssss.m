n=6
B = randn(n,n)
C = randn(n,n)
% A=zeros(n)
del=0.000000001

% B = [-1.72046870082285,-0.840848006262920,0.136059626849985,1.12596010948328,-1.49164980783487,-0.666654171832666;0.198754840450701,-0.0266104816740556,-0.542969960862893,-1.49109570625847,-1.09866949938325,0.738423128069416;1.60430658415309,0.262225204747552,0.189833529991879,1.29936401163752,0.200056167886774,0.894381149570428;-1.18456881712455,1.09948308991036,-1.76814421953103,-1.66413480558299,0.961091252946386,0.0782103389231284;1.01072120165614,1.09534226685426,0.526076301075762,1.26816370406909,-0.772756154279721,1.45081882824035;1.37739457118937,0.128533783916308,2.30392294101240,0.0556176629818286,-1.12638295737470,-0.119655803484486] 

%%%%%%%%%%%%%%%%%%% ????????????
% C=zeros(n)
% for i=-1:1
%     C=C+diag(diag(B,i),i)
% end
[Q,r]=qr(C);
d=[10,5,5.00001,3,3.001,1];
D=diag(d);
A=Q*D*Q'
A=hess(A)
% H=h_trid(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%% несемтрична€
% 
% for i=-1:1
%     A=A+diag(diag(B,i),i)
% end

%%%%%%%%% симетрична€
% A=diag(diag(B))
% A=A+diag(diag(B,1),1)
% A=A+diag(diag(B,1),-1)

if det(A)==0
    display("DET")
end
%%%%%%%%%
syms x;
ff=det(A-x*eye(n))

f=expand(dl(n,A))

R=roots(poly(A))
L_g=[]
R_g=[]
for i=1:length(R)
    R_g(i)=R(i)+0.0002
    L_g(i)=R(i)-0.0002
end
X=[]
for i=1:length(R)
   a=L_g(i);
    x0=R_g(i);
    
    x1=x0+del
    y0 = subs(f,x,x0);
    y1 = subs(f,x,x1);
    eps=10^(-4);
    n_iter=0;
%     f=subs(f,x); %функци€ численно
    if subs(f,x,a)*subs(f,x,x0)<0 %переводим функцию в сабс подставл€€ в х а
         while abs(x1-x0)>eps
           t = x1-(x1-x0)*y1/(y1-y0) %формула дл€ метода секущих
           y = subs(f,x,t)
           x0=x1; %назначаем вторую точку начальной
           y0=y1;
           x1=t; 
           y1=y;
           i=i+1;
         end

    end
     X(i)=x1  
    
end
X
if det(A)==0
    display("DET")
end
GGGG=roots(poly(A))
figure
hold on
x6 = -10:0.05:35;
axis([-0, 11, -1, 1])
plot(x6, double(subs(f, x6)))
plot(X,0, '*')

figure
hold on
x6 = -10:0.05:35;
axis([-0, 11, -1, 1])
plot(x6, double(subs(ff, x6)))
plot(GGGG,0,'o')

if det(A)==0
    display("DET")
end

function detl1=dl(n,D)
    
    syms x;
    a=diag(D,-1);
    b=diag(D);
    c=diag(D,1);
    if n==0
        detl = 1;
    elseif n==1
        detl = b(1)-x;
    else
        detl = (-1)^(2*n)*(b(n)-x)*dl(n-1,D)+(-1)^(2*n-1)*a(n-1)*(-1)^(2*n-2)*c(n-1)*dl(n-2,D);
    end
    detl1=detl;
end
function Anew = h_trid(A)
%  H_TRID(A) uses Householder method to form a tridiagonal matrix from A.
%  Must have a SQUARE SYMMETRIC matrix as the input.
%
%
%  Example:   
%
%             B=[0 1 1;1 2 1;1 1 1];   
%             h_trid(B)
%             
%
% Author:  Matt Fig
% Contact: popkenai@yahoo.com
[M N] = size(A);
if M~=N || ~isequal(A,A')  %This just screens matricies that can't work.
   error('Matrix must be square symmetric only, see help.');
end
lngth = length(A);  % Preallocations. 
v = zeros(lngth,1);  
I = eye(lngth);  
Aold = A;  
for jj=1:lngth-2  % Build each vector j and run the whole procedure.
    v(1:jj) = 0;
    S = ss(Aold,jj);
    v(jj+1) = sqrt(.5*(1+abs(Aold(jj+1,jj))/(S+2*eps)));
    v(jj+2:lngth) = Aold(jj+2:lngth,jj)*sign(Aold(jj+1,jj))...
                   /(2*v(jj+1)*S+2*eps);
    P = I-2*v*v';
    Anew = P*Aold*P;
    Aold = Anew;
end
end
function anss = ss(A,jj)
% Subfunction for h_trid.
anss = sqrt(sum(A(jj+1:end,jj).^2));
end
