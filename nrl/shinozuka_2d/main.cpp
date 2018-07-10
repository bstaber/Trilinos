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
    mesh Mesh(Comm,mesh_file,1.0);
    Epetra_Map StandardMap(-1,Mesh.n_local_nodes_without_ghosts,&Mesh.local_nodes_without_ghosts[0],0,Comm);

    //int order = Teuchos::getParameter<int>(paramList->sublist("Shinozuka"), "order");
    double L1 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "lx");
    double L2 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "ly");
    double pa = 2.0*M_PI*60.0/360.0;

    int nmc = 100000;

    for (int order=1; order<=1000; order+=10){

    Epetra_MultiVector V(StandardMap,nmc);

    for (int real=0; real<nmc; ++real){

      Teuchos::RCP<shinozuka_2d> RandomField = Teuchos::rcp(new shinozuka_2d(order,L1,L2));
      RandomField->rotation = 0.0;
      //Epetra_MultiVector V(StandardMap,1,"true");

      RandomField->rng.seed(real);
      RandomField->generator(*V(real),Mesh);

    }

    std::string path = "/home/s/staber/Trilinos_results/nrl/shinozuka_2d/";
    int NumTargetElements = 0;
    if (Comm.MyPID()==0){
        NumTargetElements = Mesh.n_nodes;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,Comm);
    Epetra_Export ExportOnRoot(StandardMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,nmc,true);

    lhs_root.PutScalar(0.0);
    lhs_root.Export(V,ExportOnRoot,Insert);
    std::string filename = path + "shinozuka_2d_layer_" + std::to_string(order) + "mtx";
    int error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);

    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
