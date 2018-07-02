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

#include "Compressible_Mooney_Transverse_Isotropic_Random_Field.hpp"

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

    Teuchos::RCP<TIMooney_RandomField> interface  = Teuchos::rcp(new TIMooney_RandomField(Comm,*paramList));

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
    Epetra_Vector mu5_gmrf(GaussMap);
    Epetra_Vector x_coord(GaussMap);
    Epetra_Vector y_coord(GaussMap);
    Epetra_Vector z_coord(GaussMap);

    int node;
    double xi, eta, zeta;
    Epetra_SerialDenseMatrix matrix_X(3,interface->Mesh->el_type);
    Epetra_SerialDenseVector vector_X(3);
    Epetra_SerialDenseVector N(interface->Mesh->el_type);
    Epetra_IntSerialDenseVector seeds(5);
    seeds(0) = 0;
    seeds(0) = 1;
    seeds(0) = 2;
    seeds(0) = 3;
    seeds(0) = 4;

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
            hexa8::shape_functions(N,xi,eta,zeta);

            vector_X.Multiply('N','N',1.0,matrix_X,N,0.0);
            interface->RandomFieldGenerator(seeds);
            interface->get_material_parameters(e_lid,gp);

            mu1_gmrf[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = interface->mu1;
            mu2_gmrf[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = interface->mu2;
            mu3_gmrf[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = interface->mu3;
            mu4_gmrf[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = interface->mu4;
            mu5_gmrf[GaussMap.LID(int(e_gid*n_gauss_cells+gp))] = interface->mu5;

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
    std::string filename = "/home/s/staber/Trilinos_results/nrl/randomfields/mu1_gmrf.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(mu2_gmrf,ExportOnRoot,Insert);
    filename = "/home/s/staber/Trilinos_results/nrl/randomfields/mu2_gmrf.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(mu3_gmrf,ExportOnRoot,Insert);
    filename = "/home/s/staber/Trilinos_results/nrl/randomfields/mu3_gmrf.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(mu4_gmrf,ExportOnRoot,Insert);
    filename = "/home/s/staber/Trilinos_results/nrl/randomfields/mu4_gmrf.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(mu5_gmrf,ExportOnRoot,Insert);
    filename = "/home/s/staber/Trilinos_results/nrl/randomfields/mu5_gmrf.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(x_coord,ExportOnRoot,Insert);
    filename = "/home/s/staber/Trilinos_results/nrl/randomfields/x_coord.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(y_coord,ExportOnRoot,Insert);
    filename = "/home/s/staber/Trilinos_results/nrl/randomfields/y_coord.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(z_coord,ExportOnRoot,Insert);
    filename = "/home/s/staber/Trilinos_results/nrl/randomfields/z_coord.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;

}
