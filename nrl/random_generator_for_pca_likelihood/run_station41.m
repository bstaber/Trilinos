%
%Brian Staber (brian.staber@gmail.com)
%

clc
clearvars
close all

modelParameters.mu   = 1e3*[1.7710, 0.0658, 0.0680, 1.4152, 0.0718];
modelParameters.beta = [25.4185, 0.0432];

optimParameters.station = 41;
optimParameters.np      = 28;

load('/home/s/staber/Trilinos_results/nrl/data/eij.mat');

angle_to_id = [5,6; 1,4; 2,3; 7,8];
Yexpi = cell(4,1);
for j = 1:4
    ID = angle_to_id(j,:);
    for k = 1:2
        Exx = [exx{ID(k)}{1}, exx{ID(k)}{2}];
        Eyy = [eyy{ID(k)}{1}, eyy{ID(k)}{2}];
        Exy = [exy{ID(k)}{1}, exy{ID(k)}{2}];
        Yexpi{j}(:,k) = log(sum(Exx.^2 + Eyy.^2 + 2*Exy.^2,1));
    end
end

fd = fopen(strcat('/home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likelihood/station',num2str(optimParameters.station),'/output.txt'),'w');
[ln,lt] = meshgrid(1e-2*(5:5:25)*sqrt(50^2+25^2), 1e-2*(4:4:20)*sqrt(50^2+25^2));
ln = ln(:);
lt = lt(:);

for k = 1:length(ln)
    modelParameters.lc    = [ln(k), lt(k)];
    modelParameters.delta = repmat(0.1,1,4);

    optimParameters.tol   = 1e-6;
    optimParameters.nmc   = 100;

    output{k} = costFunction(modelParameters,optimParameters,Yexpi);
    fprintf(fd,'%d \t %f \t %f \t %f \t %f\n',k,ln(k),lt(k),0.1,output{k}.fval);
    output{k}.ln = ln(k);
    output{k}.lt = lt(k);
    output{k}.delta = 0.1;
    output{k}.nmc = 100;
    save(strcat('result_station',num2str(optimParameters.station),'_30_07_2018.mat'),'output','-v7.3');
end
