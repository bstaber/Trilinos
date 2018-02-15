clc
clearvars
close all

[ln,lt,delta] = meshgrid(1e-2*(5:5:20)*sqrt(50^2+25^2), ...
                         1e-2*(4:4:16)*sqrt(50^2+25^2), ...
                         0.1:0.1:0.4);
                     
ln = ln(:); lt = lt(:); delta = delta(:);

modelParameters.mu    = 1e3*[1.7212, 0.0426, 0.0429, 1.3138, 0.0609];
modelParameters.beta  = [27.9525, 0.306];

output = cell(length(ln));
for k = 1:length(ln)    
    modelParameters.lc    = [ln(k), lt(k)];
    modelParameters.delta = repmat(delta(k),1,4);

    optimParameters.tol   = 1e-3;
    optimParameters.nmc   = 25;
    
    output{k} = costFunction(modelParameters,optimParameters);
    fprintf('%f \t ?f \t %f \t %f\n',ln(k),lt(k),delta(k),output{k}.fval);
end