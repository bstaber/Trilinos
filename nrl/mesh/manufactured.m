clc
clear all
close all

mesh = readgmsh('composite_hexa_man.msh');
p = mesh.POS;
t = mesh.HEXAS;

u = [0.01.*(p(:,1)-25).*p(:,2), 0.1*sin(0.5*p(:,2)).^2, cos(0.01*sign(p(:,1)-25).*(p(:,1)-25).*p(:,2)).*0.1.*p(:,2)];