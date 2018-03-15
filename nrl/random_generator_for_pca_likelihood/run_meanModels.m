clc
clearvars
close all

X = [1.847392067433885   1.603595285597212   1.889059452191638   1.705966108463277   1.901445966488945
     0.005290177619018   0.115091060136008   0.174281962380900   0.057953830846869   0.000000043004127
     0.005249557258118   0.123146429951659   0.190985328688557   0.059137183799440   0.000000043003909
     1.567800176963924   1.226476292985066   0.901021720470906   1.314269326997006   1.688723579138751
     0.079517935551905   0.054836142313959   0.035636957260321   0.061309392571596   0.110441434101412
     24.274270284627807  29.805461128950615  39.058797252543307  27.971003875518885  20.000000045959482
     0.002714534579939   0.092591836219232   0.138719816194282   0.039406836756250   0.096835629239366];

optimParameters.station = 44;
optimParameters.np      = 32;
optimParameters.tol     = 1e-6;
optimParameters.nmc     = 1;

modelParameters.lc    = [10, 5];
modelParameters.delta = repmat(delta,1,4);

output = cell(length(ln),1);

fd = fopen(strcat('/home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likelihood/station', ...
    num2str(optimParameters.station),'/output.txt'),'w');
for k = 1:size(X,2)
    modelParameters.mu   = 1e3*[X(1,k), X(2,k), X(3,k), X(4,k), X(5,k)];
    modelParameters.beta = [X(6,k), X(7,k)];
    
    output{k} = costFunction(modelParameters,optimParameters);
    fprintf(fd,'%d \t %f \t %f \t %f \t %f\n',k,ln(k),lt(k),delta(k),output{k}.fval);
    save(strcat('result_meanModels_station',num2str(optimParameters.station),'.mat'),'output','-v7.3');
end
fclose(fd);
