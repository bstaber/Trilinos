clc
clear all
close all

%% QUADRATURE

gen_laguerre_rule(nq,k1-1,0,1,'genLagRuleOutput');
x1 = dlmread('genLagRuleOutput_x.txt');
w1 = dlmread('genLagRuleOutput_w.txt');

delete('genLagRuleOutput_x.txt', 'genLagRuleOutput_w.txt');

gen_laguerre_rule(nq,k2-1,0,1,'genLagRuleOutput');
x2 = dlmread('genLagRuleOutput_x.txt');
w2 = dlmread('genLagRuleOutput_w.txt');

delete('genLagRuleOutput_x.txt', 'genLagRuleOutput_w.txt');

[x3, w3] = jacobi_rule(nq, y1, y2, 0, 1, 'jacobiRule');

[N1,N2,N3] = meshgrid(x1,x2,x3);
[W1,W2,W3] = meshgrid(w1,w2,w3);

N1 = N1(:); N2 = N2(:); N3 = N3(:);
W1 = W1(:); W2 = W2(:); W3 = W3(:);

LAMBDA = theta1*N1; MU = theta2*N2;

c1f = 100*c1;
c2f = 100*c2; alphaf = c2f/2;
lambdaf = c1f - (2/3)*c2f; 

%% SOLVER

F = [1.5, 0.3; 
     0.3, 1.5]; 
 
H = F-eye(2);

for nQ = 1:50
    
WEFF = convergenceAnalysis(nQ,H,FIXEDNODES,p,NE,GDOF,ND,...
    PHASES,sb,sc,Ae,DOFe,ISPARSE,JSPARSE,KSPARSE);

save(strcat(['WEFF_genLaguerre_n=', num2str(nQ), '.mat']), ...
    'WEFF', 'N1', 'N2', 'LAMBDA', 'MU', 'k1', 'k2', 'theta1', 'theta2', ...
    'c1', 'c2', 'c1f', 'c2f', 'lambda', 'lambdaf', 'delta1', 'delta2');

end





