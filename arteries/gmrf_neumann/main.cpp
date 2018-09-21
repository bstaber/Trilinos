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

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "neumannInnerSurface_StochasticPolyconvexHGO.hpp"
#include "newtonRaphson.hpp"
#include <boost/math/special_functions/gamma.hpp>

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

    Teuchos::RCP<neumannInnerSurface_StochasticPolyconvexHGO> my_interface =
    Teuchos::rcp(new neumannInnerSurface_StochasticPolyconvexHGO(Comm,*paramList));

    std::ifstream parameters_file_1, parameters_file_2, parameters_file_3, parameters_file_4;
    std::string path = Teuchos::getParameter<std::string>(paramList->sublist("Mesh"), "path_to_gmrf");
    parameters_file_1.open(path+"w1.txt");
    parameters_file_2.open(path+"w2.txt");
    parameters_file_3.open(path+"w3.txt");
    parameters_file_4.open(path+"w4.txt");

    unsigned int n_cells_p1_med = 297828;
    unsigned int n_nodes_p1_med = 58464;

    int error;
    std::string path_p1 = Teuchos::getParameter<std::string>(paramList->sublist("Mesh"), "path_to_p1_connectivity");

    my_interface->get_media(n_cells_p1_med,n_nodes_p1_med,path_p1);

    if (parameters_file_1.is_open() && parameters_file_2.is_open() && parameters_file_3.is_open() && parameters_file_4.is_open()){
        for (unsigned nmc=0; nmc<3; ++nmc){
            for (int i=0; i<n_nodes_p1_med; ++i){
                parameters_file_1 >> my_interface->w1_gmrf(i);
                parameters_file_2 >> my_interface->w2_gmrf(i);
                parameters_file_3 >> my_interface->w3_gmrf(i);
                parameters_file_4 >> my_interface->w4_gmrf(i);
            }
            Teuchos::RCP<newtonRaphson> Newton = Teuchos::rcp(new newtonRaphson(*my_interface,*paramList));
            Newton->setParameters(*paramList);
            Newton->Initialization();
            error = Newton->Solve_with_Aztec(true);
            if (!error){
                std::string filename1 = path + "/review/disp_realization" + std::to_string(nmc) + ".mtx";
                Newton->print_newton_solution(filename1);
                std::string filename2 = path + "/review/stress_realization" + std::to_string(nmc);
                my_interface->compute_center_cauchy_stress(*Newton->x,filename2);
            }
            else{
              if (Comm.MyPID()==0){
                std::cout << "Newton failed | realization nb. " << nmc << "\n";
              }
        }
        Comm.Barrier();
      }
        parameters_file_1.close();
        parameters_file_2.close();
        parameters_file_3.close();
        parameters_file_4.close();
  }
  else{
    std::cout << "Couldn't open one of the parameters_file.\n";
  }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
