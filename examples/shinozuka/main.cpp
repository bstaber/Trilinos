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

    //mesh_file = "/Users/Brian/Documents/Thesis/0-Trilinos/Trilinos/arteries/mesh/media_flatboundaries.msh";
    //mesh_file = "/Users/Brian/Documents/Thesis/0-Trilinos/Trilinos/nrl/mesh/composite_hexa.msh";
    
    std::string mesh_file = Teuchos::getParameter<std::string>(paramList->sublist("Shinozuka"),"mesh_file");
    mesh Mesh(Comm,mesh_file);
    Epetra_Map StandardMap(-1,Mesh.n_local_nodes_without_ghosts,&Mesh.local_nodes_without_ghosts[0],0,Comm);
    
    int order = Teuchos::getParameter<int>(paramList->sublist("Shinozuka"), "order");
    int nmc = Teuchos::getParameter<int>(paramList->sublist("Shinozuka"), "nmc");
    double L1 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "lx");
    double L2 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "ly");
    double L3 = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"), "lz");
    
    Teuchos::RCP<shinozuka> RandomField = Teuchos::rcp(new shinozuka(order,L1,L2,L3));
    //RandomField->rng.seed(std::time(0));
    
    Epetra_Vector V(StandardMap);
    //Epetra_Vector G(StandardMap);
    //Epetra_Vector B(StandardMap);
    
    double scdOrderMoment = 0.0;
    double GRFNorm2 = 0.0;
    
    if (Comm.MyPID()==0){
        std::cout << "\n E(||V||^2)/npoints = ";
    }
    for (unsigned int j=1; j<=nmc; ++j){
        RandomField->rng.seed(j);
        RandomField->generator(V,Mesh);
        V.Norm2(&GRFNorm2);
        scdOrderMoment += (double(j-1.0)/double(j))*scdOrderMoment + (1.0/double(j))*GRFNorm2;
        if (Comm.MyPID()==0){
            std::cout << scdOrderMoment/V.GlobalLength() << "\n";
        }
    }
    
    scdOrderMoment = scdOrderMoment/V.GlobalLength();
    if (Comm.MyPID()==0){
        std::cout << "\n\n E(||V||^2)/npoints = " << scdOrderMoment << "\n";
    }
    //double alpha = 1.0/(0.10*0.10); double beta = 10.0*0.10*0.10;
    //RandomField->icdf_gamma(V,G,alpha,beta);
    
    //double tau1 = 5.0; double tau2 = 8.0;
    //RandomField->icdf_beta(V,B,tau1,tau2);
    
    /*std::string path = "/Users/brian/Documents/GitHub/Trilinos_results/examples/shinozuka/";
    int NumTargetElements = 0;
    if (Comm.MyPID()==0){
        NumTargetElements = Mesh.n_nodes;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,Comm);
    Epetra_Export ExportOnRoot(StandardMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(G,ExportOnRoot,Insert);
    std::string filename = path + "shinozuka_gamma.mtx";
    int error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    
    lhs_root.PutScalar(0.0);
    lhs_root.Export(B,ExportOnRoot,Insert);
    filename = path + "shinozuka_beta.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    
    lhs_root.PutScalar(0.0);
    lhs_root.Export(V,ExportOnRoot,Insert);
    filename = path + "shinozuka_gaussian.mtx";
    error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
     */
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
