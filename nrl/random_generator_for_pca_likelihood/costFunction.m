clc
clearvars
clear all

% if (length(mu)~=5 || length(beta)~=2 || length(lc)~=2 || length(delta)~=4)
%    fprintf('Check inputs.\n'); 
% end
% writeXMLParameterList('nrl.msme.xml',mu,beta,lc,delta,nmc);

s = unix('mpirun -np 24 ./trilinos_mpi --xml-in-file="nrl.msme.xml"');
if (s~=0)
    fprintf('Trilinos program failed.\n');
    f = 0;
    return;
end

theta       = ['15';'30';'60';'75'];
theta_to_id = [5,6;1,4;2,3;7,8]; 

eta    = cell(4,1);
etaExp = cell(4,1);
m      = zeros(4,1);
for i = 1:4
    
    Y = zeros(2355,nmc);
    for j = 0:19
        path = '/Users/brian/Documents/GitHub/Trilinos_results/nrl/random_generator_for_pca_likelihood/delta=0.2_L=10percent';
        filename = strcat(path, '/RandomVariableY_angle=',num2str(theta(i,:)),'_nmc=',num2str(j),'.mtx');
        Y(:,j+1) = log(load(filename));
    end
    meanY   = mean(Y,2);
    
    for k = 1:2
        id = theta_to_id(i,k);
        filename = strcat('/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/exx_id',num2str(id),'.txt');
        exx = dlmread(filename);
        filename = strcat('/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/eyy_id',num2str(id),'.txt');
        eyy = dlmread(filename);
        filename = strcat('/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/exy_id',num2str(id),'.txt');
        exy = dlmread(filename);
        Yexp(:,k) = log(sum(exx.^2 + eyy.^2 + 2*exy.^2,1))';
    end
    
    CovarianceMatrix = cov(Y');
    [L,P]   = eig(CovarianceMatrix);
    [P,idx] = sort(diag(P),'descend');
    L       = L(:,idx);
    err     = [1; 1 - cumsum(P)/sum(P)];
    m(i)    = find(err>=1e-6,1,'last');
    
    eta{i}    = zeros(m(i),nmc);
    etaExp{i} = zeros(m(i),2);
    for l = 1:m(i)         
        eta{i}(l,:)    = ( (Y    - repmat(meanY,1,nmc))'*L(:,l) )/sqrt(P(l));
        etaExp{i}(l,:) = ( (Yexp - repmat(meanY,1,2))'  *L(:,l) )/sqrt(P(l));
    end

end

f = 0;
end
%