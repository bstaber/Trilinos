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
#include <boost/math/special_functions/gamma.hpp>

#include "shinozukapp.hpp"

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

    Teuchos::RCP<neumannInnerSurface_StochasticPolyconvexHGO> my_interface
    = Teuchos::rcp(new neumannInnerSurface_StochasticPolyconvexHGO(Comm,*paramList));

    std::ifstream parameters_file_1, parameters_file_2, parameters_file_3, parameters_file_4;

    std::string path1 = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/gmrf_neumann/a100_gamma3_delta010/";
    parameters_file_1.open(path1+"w1.txt");
    parameters_file_2.open(path1+"w2.txt");
    parameters_file_3.open(path1+"w3.txt");
    parameters_file_4.open(path1+"w4.txt");

    unsigned int n_cells_p1_med = 297828;
    unsigned int n_nodes_p1_med = 58464;

    Epetra_Map StandardMap(int(n_nodes_p1_med),0,Comm);

    Epetra_Vector c1(StandardMap);
    Epetra_Vector c2(StandardMap);
    Epetra_Vector u1(StandardMap);
    Epetra_Vector mu4(StandardMap);

    std::string path = "/Users/Brian/Documents/Thesis/Trilinos/arteries/mesh/connectivity_p1_media.txt";
    my_interface->get_media(n_cells_p1_med,n_nodes_p1_med,path);

    if (parameters_file_1.is_open() && parameters_file_2.is_open() && parameters_file_3.is_open() && parameters_file_4.is_open()){

        for (unsigned nmc=0; nmc<1; ++nmc){
            for (int i=0; i<n_nodes_p1_med; ++i){
                parameters_file_1 >> my_interface->w1_gmrf(i);
                parameters_file_2 >> my_interface->w2_gmrf(i);
                parameters_file_3 >> my_interface->w3_gmrf(i);
                parameters_file_4 >> my_interface->w4_gmrf(i);
            }
            for (int j=0; j<n_nodes_p1_med; ++j){
                if (StandardMap.MyGID(j)){
                    c1[StandardMap.LID(j)]  = my_interface->icdf_gamma(my_interface->w1_gmrf[j],my_interface->alpha1,my_interface->alpha2);
                    c2[StandardMap.LID(j)]  = my_interface->icdf_gamma(my_interface->w2_gmrf[j],my_interface->alpha3,my_interface->alpha4);
                    u1[StandardMap.LID(j)]  = my_interface->icdf_beta(my_interface->w3_gmrf[j],my_interface->tau1,my_interface->tau2);
                    mu4[StandardMap.LID(j)] = my_interface->icdf_gamma(my_interface->w4_gmrf[j],my_interface->alpha5,my_interface->alpha6);
                }
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

    int error;
    int NumTargetElements = 0;
    if (Comm.MyPID()==0){
        NumTargetElements = n_nodes_p1_med;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,Comm);
    Epetra_Export ExportOnRoot(StandardMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(c1,ExportOnRoot,Insert);
    std::string filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_nodes/c1_test.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(c2,ExportOnRoot,Insert);
    filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_nodes/c2_test.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(u1,ExportOnRoot,Insert);
    filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_nodes/u1_test.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(mu4,ExportOnRoot,Insert);
    filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_nodes/mu4_test.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);


#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
