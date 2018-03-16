clc
clearvars
close all

[ln,lt,delta] = meshgrid(1e-2*(5:5:20)*sqrt(50^2+25^2), ...
                         1e-2*(4:4:16)*sqrt(50^2+25^2), ...
                         0.1:0.1:0.4);

ln = ln(:); lt = lt(:); delta = delta(:);

modelParameters.mu      = 1e3*[1.7212, 0.0426, 0.0429, 1.3138, 0.0609];
modelParameters.beta    = [27.9525, 0.306];
optimParameters.station = 44;
optimParameters.np      = 18;

load('eij.mat');

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

output = cell(length(ln));

fd = fopen(strcat('/home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likelihood/station',num2str(optimParameters.station),'/output.txt'),'w');
for k = 1:length(ln)
    modelParameters.lc    = [ln(k), lt(k)];
    modelParameters.delta = repmat(delta(k),1,4);

    optimParameters.tol   = 1e-3;
    optimParameters.nmc   = 25;

    output{k} = costFunction(modelParameters,optimParameters,Yexpi);
    fprintf(fd,'%d \t %f \t %f \t %f \t %f\n',k,ln(k),lt(k),delta(k),output{k}.fval);
    save(strcat('result_station',num2str(optimParameters.station),'.mat'),'output','-v7.3');
end
fclose(fd);
