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
        if (Comm.MyPID()==0){
            paramList->print(std::cout,2,true,true);
        }
    }
    
    std::string mesh_file = Teuchos::getParameter<std::string>(paramList->sublist("Shinozuka"),"mesh_file");
    mesh Mesh(Comm,mesh_file);
    Epetra_Map StandardMap(-1,Mesh.n_local_nodes_without_ghosts,&Mesh.local_nodes_without_ghosts[0],0,Comm);
    
    int order = Teuchos::getParameter<int>(paramList->sublist("Shinozuka"), "order");
    double L1 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "lx");
    double L2 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "ly");
    double L3 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "lz");
    
    Teuchos::RCP<shinozuka> RandomField = Teuchos::rcp(new shinozuka(order,L1,L2,L3));
    
    Epetra_MultiVector V(StandardMap,5,"true");
    Epetra_MultiVector G(StandardMap,5,"true");
    
    RandomField->rng.seed(1);
    for (unsigned int i=0; i<5; ++i){
        RandomField->generator(*V(i),Mesh);
    }
    
    double deltaN = 0.10;
    double deltaMk = 0.10;
    
    double alpha = 3.0/(2.0*deltaN*deltaN); double beta = 1.0;
    RandomField->icdf_gamma(*V(0),*G(0),alpha,beta);
    
    alpha = 3.0/(2.0*deltaN*deltaN) - 1.0/2.0;
    RandomField->icdf_gamma(*V(1),*G(1),alpha,beta);
    
    alpha = 1.0/(deltaMk*deltaMk); beta = 1.0*deltaMk*deltaMk;
    RandomField->icdf_gamma(*V(3),*G(3),alpha,beta);
    RandomField->icdf_gamma(*V(4),*G(4),alpha,beta);
    *G(2) = *V(2);
    
    // it gives the translated field but not the coefficients M1,...,M5.
    
    std::string path = "/Users/brian/Documents/GitHub/Trilinos_results/isotransrandomfield/";
    int NumTargetElements = 0;
    if (Comm.MyPID()==0){
        NumTargetElements = Mesh.n_nodes;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,Comm);
    Epetra_Export ExportOnRoot(StandardMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,5,true);
    
    lhs_root.PutScalar(0.0);
    lhs_root.Export(G,ExportOnRoot,Insert);
    std::string filename = path + "shinozuka_translatedfield_deltaN_" + std::to_string(deltaN) + "_deltaMk_" + std::to_string(deltaMk) + ".mtx";
    int error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
     
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
