clc
clear all
close all

mesh = readgmsh('rubberblock.msh');
p = mesh.POS;
t = mesh.HEXAS;

u = load('u.mtx');
e = load('e22.mtx');
s = load('sig22.mtx');
pdef = p/1000 + [u(1:3:end), u(2:3:end), u(3:3:end)];
writevtk('rubberblock.vtk',pdef,t(:,1:8),[],[],[u(1:3:end), u(2:3:end), u(3:3:end)],[s,e]);