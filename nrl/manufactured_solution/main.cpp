#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "manufacturedSolution.hpp"

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
    
    Teuchos::RCP<manufacturedSolution> manufactured = Teuchos::rcp(new manufacturedSolution(Comm,*paramList));
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*manufactured,*paramList));
    
    int id = 1;
    Epetra_SerialDenseVector parameters(7);
    parameters(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
    parameters(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
    parameters(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
    parameters(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
    parameters(4) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
    parameters(5)  = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
    parameters(6)  = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
    for (unsigned int i=0; i<5; i++){
        parameters(i) = 1.0e9*parameters(i);
    }
    double plyagl = 2.0*M_PI*30.0/360.0;
    manufactured->set_parameters(parameters);
    manufactured->set_plyagl(plyagl);
    
    double xi = 0.0;
    Newton->Initialization();
    Newton->setParameters(*paramList);
    Newton->bc_disp = 0.0;
    double delta = 0.025;
    for (int k=1; k<5; ++k){
        manufactured->setStep(k*delta);
        int error = Newton->Solve_with_Aztec(true);
    }
    std::string path1 = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/manufactured.mtx";
    Newton->print_newton_solution(path1);
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
