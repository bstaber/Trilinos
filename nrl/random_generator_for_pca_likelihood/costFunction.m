function output = costFunction(modelParameters,optimParameters)

    if (length(modelParameters.mu)~=5 || length(modelParameters.beta)~=2 || length(modelParameters.lc)~=2 || length(modelParameters.delta)~=4)
       fprintf('Check inputs.\n');
    end
    xmlfilename = strcat('nrl.msme.station',num2str(optimParameters.station),'.xml');
    writeXMLParameterList(xmlfilename,   modelParameters.mu, ...
                                         modelParameters.beta, ...
                                         modelParameters.lc, ...
                                         modelParameters.delta, ...
                                         optimParameters.nmc, ...
                                         optimParameters.station);
    executable = strcat('mpirun -np', 32, num2str(optimParameters.np), 32, './trilinos_mpi_', ...
                         num2str(optimParameters.station), 32, '--xml-in-file="',xmlfilename,'"');
    s = unix(executable);
    if (s~=0)
        fprintf('Trilinos program failed.\n');
        return;
    end

    theta         = ['15';'30';'60';'75'];
    theta_to_id   = [5,6;1,4;2,3;7,8];

    output.eta    = cell(4,1);
    output.etaExp = cell(4,1);
    output.m      = zeros(4,1);
    output.fval   = 0;

    for i = 1:4
        Y    = zeros(2355,optimParameters.nmc);
        Yexp = zeros(2355,2);
        for j = 0:optimParameters.nmc-1
            path     = strcat('/home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likelihood/station',num2str(optimParameters.station));
            filename = strcat(path, '/RandomVariableY_angle=',num2str(theta(i,:)),'_nmc=',num2str(j),'.mtx');
            Y(:,j+1) = log(load(filename));
        end
        meanY = mean(Y,2);

        for k = 1:2
            id = theta_to_id(i,k);
            filename  = strcat('/home/s/staber/Trilinos_results/nrl/data/exx_id',num2str(id),'.txt');
            exx       = dlmread(filename);
            filename  = strcat('/home/s/staber/Trilinos_results/nrl/data/eyy_id',num2str(id),'.txt');
            eyy       = dlmread(filename);
            filename  = strcat('/home/s/staber/Trilinos_results/nrl/data/exy_id',num2str(id),'.txt');
            exy       = dlmread(filename);
            Yexp(:,k) = log(sum(exx.^2 + eyy.^2 + 2*exy.^2,1))';
        end

        CovarianceMatrix = cov(Y');
        [L,P]            = eig(CovarianceMatrix);
        [P,idx]          = sort(diag(P),'descend');
        L                = L(:,idx);
        err              = [1; 1 - cumsum(P)/sum(P)];

        output.m(i)      = find(err>=optimParameters.tol,1,'last');
        output.eta{i}    = zeros(output.m(i),optimParameters.nmc);
        output.etaExp{i} = zeros(output.m(i),2);

        for l = 1:output.m(i)
            output.eta{i}(l,:)    = ( (Y    - repmat(meanY,1,optimParameters.nmc))'*L(:,l) )/sqrt(P(l));
            output.etaExp{i}(l,:) = ( (Yexp - repmat(meanY,1,2)                  )'*L(:,l) )/sqrt(P(l));
            [~,supp]              = ksdensity(output.eta{i}(l,:));
            pts                   = find(output.etaExp{i}(l,:)>=min(supp) & output.etaExp{i}(l,:)<=max(supp));
            if (isempty(pts)==0)
                [pdf,~]           = ksdensity(output.eta{i}(l,:),output.etaExp{i}(l,pts));
                output.fval       = output.fval + sum(log(pdf));
            end
        end
    end
    unix(strcat('rm /home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likeliehood/station',num2str(optimParameters.station),'/*'));
end
