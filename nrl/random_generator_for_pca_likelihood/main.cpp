#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Newton_Raphsonpp.hpp"
#include "RandomSearch_DeterministicModel.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>

int main(int argc, char *argv[]){

    std::string xmlInFileName = "";

    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setDocString("Run: mpirun -np 28 ./trilinos_mpi --xml-in-file='nrl.aztec.linux.xml'");

    Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;
        return parse_return;
    }

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(new Teuchos::ParameterList);
    if(xmlInFileName.length()){
        Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(*paramList));
    }
    if (Comm.MyPID()==0){
        paramList->print(std::cout,2,true,true);
    }

    //for (unsigned int test=0; test<10; ++test){
      Teuchos::RCP<RandomGeneratorForPCA_andLikelihood> obj =
      Teuchos::rcp(new RandomGeneratorForPCA_andLikelihood(Comm,*paramList));

      /*Epetra_SerialDenseVector x(7);
      x(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
      x(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
      x(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
      x(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
      x(4) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
      x(5) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
      x(6) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
      for (unsigned int i=0; i<5; i++){
          x(i) = 1.0e3*x(i);
      }

      int niter   = 200;
      double tol  = 1.0e-6;

      double fval = obj->randomsearch(x,niter,tol);
    }*/

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
}