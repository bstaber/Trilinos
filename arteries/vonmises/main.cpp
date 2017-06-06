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
#include "Newton_Raphsonpp.hpp"
#include <boost/math/special_functions/gamma.hpp>

int load_solution(std::string & filesolution, Epetra_Vector * x);

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
    std::string path = Teuchos::getParameter<std::string>(paramList->sublist("Mesh"), "path_to_gmrf");
    parameters_file_1.open(path+"w1.txt");
    parameters_file_2.open(path+"w2.txt");
    parameters_file_3.open(path+"w3.txt");
    parameters_file_4.open(path+"w4.txt");
    
    unsigned int n_cells_p1_med = 297828;
    unsigned int n_nodes_p1_med = 58464;
    
    int error;
    std::string path_p1 = Teuchos::getParameter<std::string>(paramList->sublist("Mesh"), "path_to_p1_connectivity");
    interface->get_media(n_cells_p1_med,n_nodes_p1_med,path_p1);
    
    Epetra_Vector * x = new Epetra_Vector(*interface->StandardMap);
    
    if (parameters_file_1.is_open() && parameters_file_2.is_open() && parameters_file_3.is_open() && parameters_file_4.is_open()){
        
        for (unsigned nmc=0; nmc<1; ++nmc){
            for (int i=0; i<n_nodes_p1_med; ++i){
                parameters_file_1 >> interface->w1_gmrf(i);
                parameters_file_2 >> interface->w2_gmrf(i);
                parameters_file_3 >> interface->w3_gmrf(i);
                parameters_file_4 >> interface->w4_gmrf(i);
            }
            Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*interface,*paramList));
            Newton->setParameters(*paramList);
            
            std::string filesolution = path + "disp_realization" + std::to_string(nmc) + ".mtx";
            int error = load_solution(filesolution,x);
            if (!error){
                std::string filestress = path + "stress_realization" + std::to_string(nmc);
                interface->compute_gauss_vonmises(*x,filestress);
            }
            else{
                if (Comm.MyPID()==0){
                    std::cout << "Error with EpetraExt::MatrixMarketFileToVector\n";
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
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
    
}

int load_solution(std::string & filesolution, Epetra_Vector * x){
    
    int LID;
    int error = 0;
    double data;
    x->PutScalar(0.0);
    
    std::ifstream file;
    file.open(filesolution);
    if (file.is_open()){
        for (int i=0; i<x->GlobalLength(); ++i){
            file >> data;
            if ( (x->Map()).MyGID(i) ){
                LID = (x->Map()).LID(i);
                x->ReplaceMyValues(1,&data,&LID);
            }
        }
        
        std::cout << std::setw(10) << x->GlobalLength() << std::setw(10) << x->MyLength() << "\n";
        
        
    }
    else{
        std::cout << "Couldn't open the file: " << filesolution << "\n";
        error = -1;
    }
    return error;
}
