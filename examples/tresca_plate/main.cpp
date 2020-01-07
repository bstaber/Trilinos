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

#include "tresca_plate.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

int main(int argc, char *argv[]){

    std::string xmlInFileName = "";

    Teuchos::CommandLineProcessor clp(false);
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
    Teuchos::RCP<tresca_plate> problem = Teuchos::rcp(new tresca_plate(Comm,*paramList));

    std::string timestep = Teuchos::getParameter<std::string>((*paramList).sublist("Newton"), "timestep");
    if (timestep=="sequence") problem->sequence_bvp(true);
    else if (timestep=="automatic") problem->incremental_bvp(true);
    else std::cout << "Unknown timestep method (sequence, automatic)." << std::endl;

    problem->print_at_gauss_points("/Users/brian/Documents/GitHub/TrilinosUQComp/results/plasticity/plate/gausspoints.mtx");
    problem->print_mean_at_center("/Users/brian/Documents/GitHub/TrilinosUQComp/results/plasticity/plate/mean_center.mtx");

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
