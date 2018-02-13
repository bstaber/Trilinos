clc
clearvars
close all

s = unix('git pull origin master');
s = unix('make -f Makefile.linux build');