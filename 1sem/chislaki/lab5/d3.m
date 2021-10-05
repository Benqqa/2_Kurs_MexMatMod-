n=6
B = randn(n,n)
A=zeros(n)
% B = [-1.72046870082285,-0.840848006262920,0.136059626849985,1.12596010948328,-1.49164980783487,-0.666654171832666;0.198754840450701,-0.0266104816740556,-0.542969960862893,-1.49109570625847,-1.09866949938325,0.738423128069416;1.60430658415309,0.262225204747552,0.189833529991879,1.29936401163752,0.200056167886774,0.894381149570428;-1.18456881712455,1.09948308991036,-1.76814421953103,-1.66413480558299,0.961091252946386,0.0782103389231284;1.01072120165614,1.09534226685426,0.526076301075762,1.26816370406909,-0.772756154279721,1.45081882824035;1.37739457118937,0.128533783916308,2.30392294101240,0.0556176629818286,-1.12638295737470,-0.119655803484486] 

%%%%%%%%%%%%%%%%%%% ????????????
% C=zeros(n)
% for i=-1:1
%     C=C+diag(diag(B,i),i)
% end
% [Q,r]=qr(C);
% d=[10,5,5.00001,3,3.001,1];
% D=diag(d);
% A=Q.'*D*Q
%%%%%%%%%%%%%%%%%%%%%%%%%%% несемтричная
% 
% for i=-1:1
%     A=A+diag(diag(B,i),i)
% end

%%%%%%%%% симетричная
A=diag(diag(B))
A=A+diag(diag(B,1),1)
A=A+diag(diag(B,1),-1)
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
    R_g(i)=R(i)+0.002
    L_g(i)=R(i)-0.002
end
X=[]
for i=1:length(R)
    a=L_g(i);
    b=R_g(i);
    eps=10^(-4);
    n_iter=0;
    if subs(f, x,a)*subs(f, x,b)<0
        while (abs(b-a)/2)>eps
                n_iter=n_iter+1;
                c=(a+b)/2;
                if subs(f, x,a)*subs(f, x,c)<0
                    b=c;
                else
                    a=c;
                end
        end
    end
    X(i)=((a+b)/2);
    
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
dl(2,A)
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

