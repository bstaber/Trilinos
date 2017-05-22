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
    
    Teuchos::RCP<Interface_arteries> interface = Teuchos::rcp(new Interface_arteries(Comm,*paramList));
    
    std::ifstream parameters_file_1, parameters_file_2, parameters_file_3, parameters_file_4;
    
    std::string path1 = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/gmrf_neumann/a100_gamma3_delta010/";
    parameters_file_1.open(path1+"w1.txt");
    parameters_file_2.open(path1+"w2.txt");
    parameters_file_3.open(path1+"w3.txt");
    parameters_file_4.open(path1+"w4.txt");
    
    unsigned int n_cells_p1_med = 297828;
    unsigned int n_nodes_p1_med = 58464;
    
    std::string path = "/Users/Brian/Documents/Thesis/Trilinos/arteries/mesh/";
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
    
    int e_gid;
    int n_local_cells = interface->Mesh->n_local_cells;
    int n_gauss_cells = interface->Mesh->n_gauss_cells;
    std::vector<int> local_gauss_points(n_local_cells*n_gauss_cells);
    for (unsigned int e_lid=0; e_lid<n_local_cells; ++e_lid){
        e_gid = interface->Mesh->local_cells[e_lid];
        for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
            local_gauss_points[e_lid*n_gauss_cells+gp] = e_gid*n_gauss_cells+gp;
        }
        
    }
    Epetra_Map GaussMap(-1,n_local_cells*n_gauss_cells,&local_gauss_points[0],0,Comm);
    
    Epetra_Vector mu1_gmrf(GaussMap);
    Epetra_Vector mu2_gmrf(GaussMap);
    Epetra_Vector mu3_gmrf(GaussMap);
    Epetra_Vector mu4_gmrf(GaussMap);
    Epetra_Vector x_coord(GaussMap);
    Epetra_Vector y_coord(GaussMap);
    Epetra_Vector z_coord(GaussMap);

    int node;
    double xi, eta, zeta;
    Epetra_SerialDenseMatrix matrix_X(3,interface->Mesh->el_type);
    Epetra_SerialDenseVector vector_X(3);
    Epetra_SerialDenseVector N(interface->Mesh->el_type);
    for (unsigned int e_lid=0; e_lid<n_local_cells; ++e_lid){
        e_gid = interface->Mesh->local_cells[e_lid];
        for (int inode=0; inode<interface->Mesh->el_type; ++inode){
            node = interface->Mesh->cells_nodes[interface->Mesh->el_type*e_gid+inode];
            matrix_X(0,inode) = interface->Mesh->nodes_coord[3*node+0];
            matrix_X(1,inode) = interface->Mesh->nodes_coord[3*node+1];
            matrix_X(2,inode) = interface->Mesh->nodes_coord[3*node+2];
        }
        for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
            xi   = interface->Mesh->xi_cells[gp];
            eta  = interface->Mesh->eta_cells[gp];
            zeta = interface->Mesh->zeta_cells[gp];
            tetra10::shape_functions(N,xi,eta,zeta);
            vector_X.Multiply('N','N',1.0,matrix_X,N,0.0);
            interface->get_material_parameters(e_lid,gp);
            mu1_gmrf[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = interface->mu1;
            mu2_gmrf[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = interface->mu2;
            mu3_gmrf[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = interface->mu3;
            mu4_gmrf[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = interface->mu4;
            x_coord[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = vector_X(0);
            y_coord[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = vector_X(1);
            z_coord[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = vector_X(2);
         }
    }
    
    int error;
    int NumTargetElements = 0;
    if (Comm.MyPID()==0){
        NumTargetElements = interface->Mesh->n_cells*n_gauss_cells;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,Comm);
    Epetra_Export ExportOnRoot(GaussMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(mu1_gmrf,ExportOnRoot,Insert);
    std::string filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_gauss/mu1_gmrf.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(mu2_gmrf,ExportOnRoot,Insert);
    filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_gauss/mu2_gmrf.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    
    lhs_root.PutScalar(0.0);
    lhs_root.Export(mu3_gmrf,ExportOnRoot,Insert);
    filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_gauss/mu3_gmrf.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    
    lhs_root.PutScalar(0.0);
    lhs_root.Export(mu4_gmrf,ExportOnRoot,Insert);
    filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_gauss/mu4_gmrf.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    
    lhs_root.PutScalar(0.0);
    lhs_root.Export(x_coord,ExportOnRoot,Insert);
    filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_gauss/x_coord.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(y_coord,ExportOnRoot,Insert);
    filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_gauss/y_coord.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    
    lhs_root.PutScalar(0.0);
    lhs_root.Export(z_coord,ExportOnRoot,Insert);
    filename = "/Users/Brian/Documents/Thesis/Trilinos_results/arteries/boost_gmrf/gmrf_gauss/z_coord.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
    
}
