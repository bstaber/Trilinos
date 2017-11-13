#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "rubberblock.hpp"
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
    
    Epetra_SerialDenseVector parameters(7);
    Teuchos::RCP<rubberblock> interface = Teuchos::rcp(new rubberblock(Comm,*paramList));
    
    parameters(0)    = Teuchos::getParameter<double>(paramList->sublist("rubberblock"),"lambda");
    parameters(1)    = Teuchos::getParameter<double>(paramList->sublist("rubberblock"),"mu");
    interface->set_parameters(parameters);
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*interface,*paramList));
    
    if (Comm.MyPID()==0){
        std::cout << "Top = " << interface->topcoord << "\n";
    }
    
    double xi = 0.0;
    double g  = -5.0*0.3;
    std::string pathsig22 = "/Users/brian/Documents/GitHub/Trilinos/examples/rubberblock/sig22.mtx";
    std::string pathsolut = "/Users/brian/Documents/GitHub/Trilinos/examples/rubberblock/u.mtx";
    
    Newton->Initialization();
    Newton->setParameters(*paramList);
    Newton->bc_disp = g;
    int error = Newton->Solve_with_Aztec(true);
    Newton->print_newton_solution(pathsolut);
    interface->compute_cauchy(*Newton->x,xi,xi,xi,pathsig22);
    
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
