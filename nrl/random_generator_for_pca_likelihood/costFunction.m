function f = main(mu,beta,lc,delta)

if (length(mu)~=5 || length(beta)~=2 || length(lc)~=2 || length(delta)~=4)
   fprintf('Check inputs.\n'); 
end

writeXMLParameterList('nrl.msme.xml',mu,beta,lc,delta);

f = 0;
end
%s = unix('mpirun -np 24 ./trilinos_mpi --xml-in-file="nrl.msme.xml"');