clc
clearvars
close all

X = [1.847305797816567   1.604001927458218   1.907043912388506   1.706486189075632   1.769017508886223
    0.005376710535645   0.114625779806038   0.154036583534123   0.057388351625262   0.000008510008773
    0.005139929801197   0.122739679192122   0.167208494913353   0.059077990032629   0.000008510039000
    1.567795920644496   1.226504600774701   0.902150886924485   1.314292350043600   2.019537222717919
    0.079518524907316   0.054865622201899   0.036767744476560   0.061332096908722   0.147713367550041
   24.274381300141435  29.806819794191611  39.150056931099662  27.971878357799234  16.420347459340434
   0.002730236643057   0.093254046036705   0.184311975895317   0.039848773492742   0.000003791513270];

% X= [1.0;
%     1.0;
%     1.0;
%     3.479232101380515
%     1.029493689764523
%     5.042703469681112
%     0.000000000142884];

optimParameters.station = 44;
optimParameters.np      = 32;
optimParameters.tol     = 1e-6;
optimParameters.nmc     = 25;

modelParameters.lc    = [11.180339887498949, 4.472135954999580];
modelParameters.delta = repmat(0.1,1,4);

output = cell(size(X,2),1);

%fd = fopen(strcat('/home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likelihood/station', ...
%    num2str(optimParameters.station),'/output.txt'),'w');

load('eij.mat');
angle_to_id = [5,6; 1,4; 2,3; 7,8];
Yexpi = cell(4,1);
for j = 1:4
    ID = angle_to_id(j,:);
    for k = 1:2
        Exx = [exx{ID(k)}{1}, exx{ID(k)}{2}];
        Eyy = [eyy{ID(k)}{1}, eyy{ID(k)}{2}];
        Exy = [exy{ID(k)}{1}, exy{ID(k)}{2}];
        Yexp{j}(:,k) = log(sum(Exx.^2 + Eyy.^2 + 2*Exy.^2,1));
    end
end

for k = 1:size(X,2)
    modelParameters.mu   = 1e3*[X(1,k), X(2,k), X(3,k), X(4,k), X(5,k)];
    modelParameters.beta = [X(6,k), X(7,k)];

    output{k} = costFunction(modelParameters,optimParameters,Yexpi);
%    fprintf(fd,'%d \t %f \t %f \t %f \t %f\n',k,ln(k),lt(k),delta(k),output{k}.fval);
    save(strcat('result_meanModels_lnltdelta14_station',num2str(optimParameters.station),'.mat'),'output','-v7.3');
end
%fclose(fd);
