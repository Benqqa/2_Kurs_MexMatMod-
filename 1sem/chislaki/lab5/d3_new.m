n=6
% B = randn(n,n)
B = [-1.72046870082285,-0.840848006262920,0.136059626849985,1.12596010948328,-1.49164980783487,-0.666654171832666;0.198754840450701,-0.0266104816740556,-0.542969960862893,-1.49109570625847,-1.09866949938325,0.738423128069416;1.60430658415309,0.262225204747552,0.189833529991879,1.29936401163752,0.200056167886774,0.894381149570428;-1.18456881712455,1.09948308991036,-1.76814421953103,-1.66413480558299,0.961091252946386,0.0782103389231284;1.01072120165614,1.09534226685426,0.526076301075762,1.26816370406909,-0.772756154279721,1.45081882824035;1.37739457118937,0.128533783916308,2.30392294101240,0.0556176629818286,-1.12638295737470,-0.119655803484486] 

A=zeros(n)
del=0.0001

%%%%%%%%%%%%%%%%%%% ��� ����� ��������� ������ ������� �������
C=zeros(n)
for i=-1:1
C=C+diag(diag(B,i),i)
end
[Q,r]=qr(C);
d=[10,5,5.00001,3,3.001,1];
D=diag(d);
A=Q'*D*Q
%%%%%%%
% for i=-1:1 %������� ����������������
%     A=A+diag(diag(B,i),i)
% end
%%%%%%%%% ������������
% A=diag(diag(B))
% A=A+diag(diag(B,1),1)
% A=A+diag(diag(B,1),-1)
%%%%%%%%%
syms x;
% f=det(A-x*eye(n)) 
% figure
% ezplot(f)

f=dl(n,A) 

Roots=roots(poly(A))
left=[]; right=[]
for i=1:length(Roots)
    right(i)=Roots(i)+0.0002
    left(i)=Roots(i)-0.0002
end
X=[]

% 

for i=1:length(Roots)
    
    a=left(i);
    x0=right(i);
    
    x1=x0+del
    y0 = subs(f,x,x0);
    y1 = subs(f,x,x1);
    eps=10^(-4);
    n_iter=0;
%     f=subs(f,x); %������� ��������
    if subs(f,x,a)*subs(f,x,x0)<0 %��������� ������� � ���� ���������� � � �
         while abs(x1-x0)>eps
           x8 = x1-(x1-x0)*y1/(y1-y0) %������� ��� ������ �������
           y = subs(f,x,x8)
           x0=x1; %��������� ������ ����� ���������
           y0=y1;
           x1=x8; 
           y1=y;
           i=i+1;
         end

    end
     X(i)=x1  
end

X
figure
hold on
x = -6:0.05:15;
axis([-5,14,-5,5])
plot(x,double(subs(f,x)))
plot(X,0,'b*')
title('������������������ ���������')

function detl=dl(n,D)
    syms x;
    a=diag(D,-1); %������� ����������
    b=diag(D);
    c=diag(D,1);
    if n==0
        detl = 1;
    elseif n==1
        detl = b(1)-x;
    else
        detl = (-1)^(2*n)*(b(n)-x)*dl(n-1,D)+(-1)^(2*n-1)*a(n-1)*(-1)^(2*n-2)*c(n-1)*dl(n-2,D);
    end
end
