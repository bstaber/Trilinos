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
#include "NeumannInnerSurface_PolyconvexHGO.hpp"
#include "Newton_Raphsonpp.hpp"

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
    }

    if (Comm.MyPID()==0){
        paramList->print(std::cout,2,true,true);
    }

    Teuchos::RCP<NeumannInnerSurface_PolyconvexHGO> my_interface = Teuchos::rcp(new NeumannInnerSurface_PolyconvexHGO(Comm,*paramList));
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*my_interface,*paramList));

    Newton->Initialization();
    Newton->setParameters(*paramList);
    int error = Newton->Solve_with_Aztec(true);
    if (!error){
        std::string path = "/home/s/staber/Trilinos_results/arteries/gmrf_neumann/";
        std::string filename1 = path + "disp_mean_model.mtx";
        Newton->print_newton_solution(filename1);
        std::string filename2 = path + "stress_mean_model";
        my_interface->compute_center_cauchy_stress(*Newton->x,filename2);
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
