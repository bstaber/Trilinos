#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <fstream>
#include "meshpp.hpp"
#include "laplacepp.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Arteries_ModelC_deterministic_neumann.hpp"
#include "Newton_Raphsonpp.hpp"

int main(int argc, char *argv[]){
    
    std::string    xmlInFileName = "";
    
    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setDocString(
                     "This program only takes one parameter list from an XMF file"
                     " given by --xml-in-file=xmlInFileName."
                     " The parameter file should contain:"
                     " - a sublist Mesh (mandatory) with the parameters mesh_file, boundary_file, nb_phys_groups."
                     " - a sublist Laplace (optional)."
                     " - a sublist Linear Solver Builder (optional)."
                     " - a sublist Newton (mandatory)."
                     );
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
        
    Teuchos::RCP<NeumannInnerSurface_PolyconvexHGO> my_interface = Teuchos::rcp(new NeumannInnerSurface_PolyconvexHGO(Comm,*paramList));
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*my_interface,*paramList));
    Newton->Initialization();
    
    Teuchos::RCP<Teuchos::ParameterList> solverBuilderSL = Teuchos::sublist(paramList,"Linear Solver Builder",true);
    int error = Newton->Solve_with_Stratimikos(solverBuilderSL);

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    
    return 0;
}
