#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "compressibleHyperelasticity_linearPatchTest.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

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
    if(xmlInFileName.length()) {
        Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(*paramList));
    }

    if (Comm.MyPID()==0){
        paramList->print(std::cout,2,true,true);
    }

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
    double angle = 15.0*2.0*M_PI/360.0;

    Teuchos::RCP<compressibleHyperelasticity_linearPatchTest> obj =
    Teuchos::rcp(new compressibleHyperelasticity_linearPatchTest(Comm,*paramList));
    obj->set_parameters(x,angle);

    Teuchos::RCP<newtonRaphson> newton =
    Teuchos::rcp(new newtonRaphson(*obj,*paramList));

    newton->Initialization();
    int flagNewton = newton->Solve_with_Aztec(true);

    double errorL2 = 1.0;
    if (!flagNewton){
      errorL2 = obj->errorL2(*newton->x);
      std::string filenameU = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/UlinearPatchTest.mtx";
      newton->print_newton_solution(filenameU);
      if(Comm.MyPID()==0){
        std::cout << "-----------------------------\n";
        std::cout << "||u||_L^2 = " << errorL2 << "\n";
        std::cout << "-----------------------------\n";
      }
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
