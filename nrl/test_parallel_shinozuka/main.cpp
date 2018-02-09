#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "shinozukapp_2d.hpp"

int main(int argc, char *argv[]){

    std::string    xmlInFileName = "";
    std::string    extraXmlFile = "";
    std::string    xmlOutFileName = "paramList.out";

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
        if (Comm.MyPID()==0){
            paramList->print(std::cout,2,true,true);
        }
    }

    std::string mesh_file = Teuchos::getParameter<std::string>(paramList->sublist("Shinozuka"),"mesh_file");
    int order = Teuchos::getParameter<int>(paramList->sublist("Shinozuka"), "order");
    double L1 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "lx");
    double L2 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "ly");

    mesh Mesh(Comm,mesh_file);
    
    Teuchos::RCP<shinozuka_layeredcomp_2d> Generator_Shinozuka = Teuchos::rcp(new shinozuka_layeredcomp_2d(order,L1,L2));

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
