clc
clear all
close all

addpath Laguerre/
addpath Jacobi/
addpath Laguerre/polynomials/

%% QUADRATURE
order = 10;

nq = 10;
k1 = 1/100; theta1 = 100;
k2 = 1/100; theta2 = 100;
y1 = 5; y2 = 5;

GenLaguerreRule(nq,k1-1,0,1,'genLagRuleOutput');
x1 = dlmread('genLagRuleOutput_x.txt');
w1 = dlmread('genLagRuleOutput_w.txt');
delete('genLagRuleOutput_x.txt', 'genLagRuleOutput_w.txt', 'genLagRuleOutput_r.txt');

GenLaguerreRule(nq,k2-1,0,1,'genLagRuleOutput');
x2 = dlmread('genLagRuleOutput_x.txt');
w2 = dlmread('genLagRuleOutput_w.txt');
delete('genLagRuleOutput_x.txt', 'genLagRuleOutput_w.txt', 'genLagRuleOutput_r.txt');

[x3, w3] = JacobiRule(nq,y1,y2,0,1,'jacobiRule');
delete('jacobiRule_x.txt', 'jacobiRule_w.txt', 'jacobiRule_r.txt');

[N1,N2,N3] = meshgrid(x1,x2,x3);
[W1,W2,W3] = meshgrid(w1,w2,w3);

N1 = N1(:); N2 = N2(:); N3 = N3(:);
W1 = W1(:); W2 = W2(:); W3 = W3(:);

X1 = theta1*N1; 
X2 = theta2*N2;
X3 = 2*X2.*N3+X1;

phi1 = repmat(sqrt(gamma(k1)*factorial(0:order)./gamma(k1+(0:order))),length(x1),1).*lf_function(length(x1),order,k1-1,x1);
phi2 = repmat(sqrt(gamma(k1)*factorial(0:order)./gamma(k1+(0:order))),length(x2),1).*lf_function(length(x2),order,k2-1,x2);
%phi3 = ...

