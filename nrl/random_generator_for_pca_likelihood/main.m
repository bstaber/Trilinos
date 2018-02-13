clc
clearvars
close all

s = unix('mpirun -np 24 ./trilinos_mpi --xml-in-file="nrl.msme.xml');