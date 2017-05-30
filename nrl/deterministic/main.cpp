#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "NRL_ModelF.hpp"
#include "Newton_Raphsonpp.hpp"

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
    
    Teuchos::RCP<NRL_ModelF> interface = Teuchos::rcp(new NRL_ModelF(Comm,*paramList));
    double m1     = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"m1");
    double m2     = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"m2");
    double beta3  = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"beta3");
    double beta4  = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"beta4");
    double beta5  = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"beta5");
    double plyagl = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"angle");
    interface->set_parameters(m1,m2,beta3,beta4,beta5);
    interface->set_plyagl(plyagl);
    
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*interface,*paramList));
    Newton->Initialization();
    
    int error = Newton->Solve_with_Aztec();
    std::string name = "/Users/Brian/Documents/Thesis/Trilinos_results/nrl/deterministic/Newton_solution.mtx";
    Newton->print_newton_solution(name);
    
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
