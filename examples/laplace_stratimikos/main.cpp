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

int main(int argc, char *argv[]){
    
    std::string    xmlInFileName = "";
    std::string    extraXmlFile = "";
    std::string    xmlOutFileName = "paramList.out";
    
    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setOption("xml-out-file",&xmlOutFileName,"The XML file to write the final parameter list to");
    clp.setDocString(
                     "This example program shows how to read in a parameter list from an"
                     " XML file (given by --xml-in-file=xmlInFileName)."
                     " The final parameter list is then written back to an XML file."
                     " (given by --xml-out-file=xmlOutFileName)."
                     "This program also shows how to use Stratikimos for a"
                     " sequence of two linear systems."
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
    
    Teuchos::RCP<Teuchos::ParameterList> solverBuilderSL = Teuchos::sublist(paramList,"Linear Solver Builder",true);
    
    Teuchos::RCP<laplace> Laplace = Teuchos::rcp(new laplace(Comm,paramList->sublist("Mesh")));
    
    Teuchos::RCP<Epetra_FECrsMatrix> matrix = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*Laplace->FEGraph));
    Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(*Laplace->StandardMap));
    Teuchos::RCP<Epetra_FEVector> rhs = Teuchos::rcp(new Epetra_FEVector(*Laplace->StandardMap));
    Teuchos::RCP<Epetra_MultiVector> rhsv = rhs;
    
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(solverBuilderSL);
    
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
    linearSolverBuilder.createLinearSolveStrategy("");
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows = lowsFactory->createOp();
    
    Thyra::SolveStatus<double> status;
    
    Teuchos::RCP<const Thyra::LinearOpBase<double>> A = Thyra::epetraLinearOp(matrix);
    Teuchos::RCP<Thyra::VectorBase<double>> x = Thyra::create_Vector(lhs,A->domain());
    Teuchos::RCP<const Thyra::MultiVectorBase<double>> b = Thyra::create_MultiVector(rhsv,A->range());
        
    int bc_indx[2];
    double bc_val[2];
    bc_indx[0] = 2; bc_indx[1] = 3;
    bc_val[0] = 0.0; bc_val[1] = 1.0;
    
    for (unsigned iter=0; iter<5; ++iter){
    Laplace->assembling_OAZ(*matrix, *rhs, &bc_indx[0], &bc_val[0]);
    rhsv = rhs;
    lhs->PutScalar(0.0);
    
    Thyra::initializeOp<double>(*lowsFactory, A, lows.ptr());
    status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *b, x.ptr());
    if (status.solveStatus){
        if (Comm.MyPID()==0){
            std::cout << "STRATIMIKOS FAILED TO CONVERGE.\n";
            std::cout << status.message << "\n";
        }
    }
    std::string name = "belos_test_" + std::to_string(iter) + ".mtx";
    int error = Laplace->print_solution(*lhs,name);
    }
    
    /*//if minor changes and reuse factorization information
    Thyra::uninitializeOp<double>(*lowsFactory, lows.ptr());
    Laplace->assembling_OAZ(*matrix, *rhs, &bc_indx[0], &bc_val[0]);
    rhsv = rhs;
    lhs->PutScalar(0.0);
    
    Thyra::initializeAndReuseOp<double>(*lowsFactory, A, lows.ptr());
    status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *b, x.ptr());
    if (status.solveStatus){
        if (Comm.MyPID()==0){
            std::cout << "STRATIMIKOS FAILED TO CONVERGE.\n";
            std::cout << status.message << "\n";
        }
    }
    error = Laplace->print_solution(*lhs,"belos_test_2.mtx");
    
    //if major changes
    bc_indx[0] = 0; bc_indx[1] = 1;
    Thyra::uninitializeOp<double>(*lowsFactory, lows.ptr());
    Laplace->assembling_OAZ(*matrix, *rhs, &bc_indx[0], &bc_val[0]);
    rhsv = rhs;
    Thyra::initializeOp<double>(*lowsFactory, A, lows.ptr());
    lhs->PutScalar(0.0);
    status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *b, x.ptr());
    
    if (status.solveStatus){
        if (Comm.MyPID()==0){
            std::cout << "STRATIMIKOS FAILED TO CONVERGE.\n";
            std::cout << status.message << "\n";
        }
    }
    error = Laplace->print_solution(*lhs,"belos_test_3.mtx");*/
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    
    return 0;
}
