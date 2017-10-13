#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "NRL_ModelF.hpp"
#include "Newton_Raphsonpp.hpp"
#include "objectiveFunction.hpp"
#include <random>

int main(int argc, char *argv[]){
    
    std::string    xmlInFileName = "";
    std::string    xmlExtraFileName = "";
    
    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setOption("extra-file",&xmlExtraFileName,"Extra XML file to read into a parameter list");
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
    if(xmlExtraFileName.length()) {
        Teuchos::updateParametersFromXmlFile(xmlExtraFileName, inoutArg(*paramList));
    }
    if (Comm.MyPID()==0){
        paramList->print(std::cout,2,true,true);
    }
    
    Teuchos::RCP<objectiveFunction<double>> obj = Teuchos::rcp(new objectiveFunction<double>(Comm,*paramList));
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
