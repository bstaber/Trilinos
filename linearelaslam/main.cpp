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
    
    rhs.PutScalar(0.0);
    lhs.PutScalar(0.0);
    
    interface->assemble_dirichlet(stiffness);
    interface->apply_dirichlet_conditions(stiffness,rhs,displacement);
    
    Epetra_LinearProblem problem;
    AztecOO solver;
    
    problem.SetOperator(&stiffness);
    problem.SetLHS(&lhs);
    problem.SetRHS(&rhs);
    solver.SetProblem(problem);
    solver.SetParameters(*paramList->sublist("Krylov"));
    
    //solver.Iterate(2000,1e-6);
    
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
