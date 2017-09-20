#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_RCP.hpp"
#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosEpetraAdapter.hpp"
#include <BelosSolverFactory.hpp>
#include "BelosBlockGmresSolMgr.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "orthotropicRF_laminate.hpp"

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
    
    Teuchos::RCP<OrthotropicRF_Laminate> interface = Teuchos::rcp(new OrthotropicRF_Laminate(Comm,*paramList));
    
    double displacement = 1.0/1000.0;
    interface->dead_pressure(0) = 0.0;
    interface->dead_pressure(1) = 1.0;
    interface->dead_pressure(2) = 0.0;
    
    Epetra_FECrsMatrix stiffness(Copy,*interface->FEGraph);
    Epetra_FEVector rhs(*interface->StandardMap);
    Epetra_Vector lhs(*interface->StandardMap);
    
    interface->material_stiffness_and_rhs_dirichlet(stiffness);
    
    //interface->assemble_dirichlet_dead_neumann(stiffness,rhs);
    //interface->apply_dirichlet_conditions(stiffness,rhs,displacement);
    
    //->setup_bcs(choose displacement)
    //get_lhs_and_rhs
    //apply_bcs
    //solve
    
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
