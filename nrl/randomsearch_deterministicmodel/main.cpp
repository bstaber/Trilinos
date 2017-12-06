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

    std::string    xmlInFileName = "";

    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setDocString("TO DO.");

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

    Teuchos::RCP<RandomSearch_DeterministicModel> obj = Teuchos::rcp(new RandomSearch_DeterministicModel(Comm,*paramList));
    
    Epetra_SerialDenseVector x(7);
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
    
    int id = 0;
    int niter = 1000;
    std::cout << "DEBUG\n";
    double fval = obj->randomsearch(x,id,niter);
    x = obj->solution;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
}

