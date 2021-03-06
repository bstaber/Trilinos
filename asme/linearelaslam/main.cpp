/*
Brian Staber (brian.staber@gmail.com)
*/

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
#include "asmeSBVP.hpp"

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

    Teuchos::RCP<asmeSBVP> interface = Teuchos::rcp(new asmeSBVP(Comm,*paramList));

    double deltaN  = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "deltaN");
    double deltaM4 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "deltaM4");
    double deltaM5 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "deltaM5");
    std::string path = Teuchos::getParameter<std::string>(paramList->sublist("Mesh"), "path");

    interface->_deltaN  = deltaN;
    interface->_deltaM4 = deltaM4;
    interface->_deltaM5 = deltaM5;

    int * seed = new int [5];
    double displacement = 1.0;

    for (unsigned int j=0; j<5000; ++j){
        seed[0] = 5*j; seed[1] = 5*j+1; seed[2] = 5*j+2; seed[3] = 5*j+3; seed[4] = 5*j+4;
        interface->solveOneRealization(displacement,seed);
        std::string path1 = path + std::to_string(j) + ".mtx";
        std::string path2 = path + std::to_string(j);
        interface->print_solution(path1);
        interface->recover_cauchy_stress(path2,seed);
        interface->compute_deformation(*interface->solution,path2,true,false);
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
