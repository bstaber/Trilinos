clc
clear all
close all

mesh = readgmsh('composite_hexa_man.msh');
p = mesh.POS;
t = mesh.HEXAS;

u = [1e-4.*(p(:,1)-25).*p(:,2).*(25-p(:,2)), 1e-3*p(:,2).*(25-p(:,2)), sin(1e-3.*(p(:,1)-25).*p(:,2).*(25-p(:,2)))];

syms a b real
syms mu1 mu2 mu3 mu4 mu5 d beta4 beta5 real

x = sym('x%d',[3,1]);
u = sym('u%d',[3,1]);
assume(x,'real');
assume(u,'real');
u(1) = a*(x(1)-25)*x(2)*(25-x(2));
u(2) = a*x(2)*(25-x(2));
u(3) = ((b/a)*u(1));

for i = 1 : 3
    for j = 1 : 3
        gradu(i,j) = diff(u(i),x(j));
    end
end
F = eye(3)+gradu;
C = F'*F;

syms c s real
n = [c;s;0];
P = n*n';
Q = [s*s, -c*s, 0; -c*s, c*c, 0; 0, 0, 1];
M = mu4*P + mu5*Q;

I1 = trace(C);
I2 = 0.5*(I1^2-trace(C*C));
I3 = det(C);
J4 = trace(C*M);
J5 = trace(I3*inv(C)*M);
g = trace(M);
Siso = 2*mu1*eye(3) + 2*mu2*(I1*eye(3)-C) + (2*mu3*I3 - d)*inv(C);
Sani = (2*J4^beta4/g^beta4)*M + (2*J5^beta5/g^beta5)*(J5*inv(C)-I3*inv(C)*M*inv(C)) - 2*g*I3^(-1/2)*inv(C);
S = Siso + Sani;
P = F*S;

f(1) = - diff(P(1,1),x(1)) - diff(P(1,2),x(2)) - diff(P(1,3),x(3));
f(2) = - diff(P(2,1),x(1)) - diff(P(2,2),x(2)) - diff(P(2,3),x(3));
f(3) = - diff(P(3,1),x(1)) - diff(P(3,2),x(2)) - diff(P(3,3),x(3));