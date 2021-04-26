f_ii= @(x,y,z)2*x.*z+2*y-4*x;
true_f=@(x) x + exp(x.^2);
a = 0;
b = 1;
n = 4;
h =(b-a)/(n);
p = @(x) -2*x;
q = @(x) -2;
f = @(x) -4*x;

y = zeros(1,n);

y0a = 1;
dy0a = 1;
y0b = 1 + exp(1);
dy0b = 1 + 2*exp(1);

alfa_0 = true_f(a);

alfa_1 = 3;
alfa_1 = 0;

beta_0 = true_f(b);

beta_1 = 1;
beta_1 = 0;
A = alfa_0*y0a + alfa_1*dy0a
B = beta_0*y0b + beta_1*dy0b
x = []
for i = 1:n
  x(i) = a + h*(i - 1);
end

ai3 = @(x)  (1 + p(x)*h/2);
ai2 = @(x)  ((h^2)*q(x) - 2);
ai1 = @(x)  (1-p(x)*h/2);
bi = @(x)  (h^2)*f(x);
% alfa_0*y(0)+alfa_1*(y(1)-y(0))./h=A
% y(0)*(alfa_0-alfa_1./h)+alfa_1*y(1)./h=A
%y(2)=(A-y0a*(alfa_0-alfa_1./h))./(h*alfa_1)

P(1) = -alfa_1./(h*(alfa_0-(alfa_1./h)))
Q(1) = A./(alfa_0-(alfa_1./h))

P(1) = -0/alfa_0;
Q(1) =alfa_0/alfa_0

for i=2:(n-1)
    P(i)= -ai3(x(i+1))/(ai1(x(i-1))*P(i-1)+ai2(x(i)))
    Q(i)=(bi(x(i))-ai1(x(i-1))*Q(i-1))/(ai1(x(i-1))*P(i-1)+ai2(x(i)))
end

% P(n-1) = (beta_0+(beta_1./h)).*h./beta_1
% Q(n-1) = -B.*h/(beta_1)
% y(n)=(bi(x(n))-ai1(x(n-1))*Q(n-1))/(ai2(x(n))+ai1(x(n-1))*P(n-1))
% y(n)=(B+(beta_1./h)*Q(n-1))/(-beta_1./h*Q(n-1)+(beta_0+(beta_1./h)))

y(n)=beta_0/alfa_0

for i=(n-1):-1:1;
    y(i)=(P(i)*y(i+1))+Q(i)
end


figure
hold on
grid on
plot(x,y,'-*b','Linewidth',1.2)
plot(x,true_f(x),'-*r','Linewidth',1.2)
