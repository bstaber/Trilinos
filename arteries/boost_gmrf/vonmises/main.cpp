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
#include "Arteries_ModelC_gmrf_neumann.hpp"

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
    
    Teuchos::RCP<Interface_arteries> interface = Teuchos::rcp(new Interface_arteries(Comm,*paramList));
    
    std::ifstream parameters_file_1, parameters_file_2, parameters_file_3, parameters_file_4;
    
    std::string path1 = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/gmrf_neumann/a100_gamma3_delta010/";
    parameters_file_1.open(path1+"w1.txt");
    parameters_file_2.open(path1+"w2.txt");
    parameters_file_3.open(path1+"w3.txt");
    parameters_file_4.open(path1+"w4.txt");
    
    unsigned int n_cells_p1_med = 297828;
    unsigned int n_nodes_p1_med = 58464;
    
    std::string path = "/Users/Brian/Documents/Thesis/Trilinos/arteries/mesh/connectivity_p1_media.txt";
    interface->get_media(n_cells_p1_med,n_nodes_p1_med,path);
    
    if (parameters_file_1.is_open() && parameters_file_2.is_open() && parameters_file_3.is_open() && parameters_file_4.is_open()){
        
        for (unsigned nmc=0; nmc<1; ++nmc){
            for (int i=0; i<n_nodes_p1_med; ++i){
                parameters_file_1 >> interface->w1_gmrf(i);
                parameters_file_2 >> interface->w2_gmrf(i);
                parameters_file_3 >> interface->w3_gmrf(i);
                parameters_file_4 >> interface->w4_gmrf(i);
            }
        }
        Comm.Barrier();
        parameters_file_1.close();
        parameters_file_2.close();
        parameters_file_3.close();
        parameters_file_4.close();
    }
    else{
        std::cout << "Couldn't open one of the parameters_file.\n";
    }
    
    Epetra_Vector x(*interface->StandardMap);
    std::ifstream displacement;
    std::string filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/gmrf_neumann/a100_gamma3_delta010/disp_realization0.mtx";
    displacement.open(filename.c_str());
    double input;
    if (displacement.is_open()){
        for (int i=0; i<3*interface->Mesh->n_nodes; ++i){
            displacement >> input;
            if (interface->StandardMap->MyGID(i)){
                x[interface->StandardMap->LID(i)] = input;
            }
        }
    }
    else{
        std::cout << "Couldn't open the displacement file.\n";
    }
    
    std::string filename_out = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/vonmises/stress";
    interface->compute_mean_cauchy_stress(x,filename_out);
    
    filename_out = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/vonmises/disp.mtx";
    int NumTargetElements = 0;
    if (Comm.MyPID()==0){
        NumTargetElements = 3*interface->Mesh->n_nodes;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,Comm);
    Epetra_Export ExportOnRoot(*interface->StandardMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(x,ExportOnRoot,Insert);
    
    int error = EpetraExt::MultiVectorToMatrixMarketFile(filename_out.c_str(),lhs_root,0,0,false);
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
    
}
