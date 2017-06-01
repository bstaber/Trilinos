#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "NRL_RandomField.hpp"
#include "Newton_Raphsonpp.hpp"

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
    
    Epetra_IntSerialDenseVector seeds(2); seeds(0) = 0; seeds(1) = 1;
    Epetra_SerialDenseVector mean(2); mean(0) = 5000.0; mean(1) = 7000.0;
    Epetra_SerialDenseVector delta(2); delta(0) = 0.4; delta(1) = 0.4;
    Epetra_SerialDenseVector corrlength(3); corrlength(0) = (50.0/1000.0)*0.50; corrlength(1) = (25.0/1000.0)*0.50; corrlength(2) = 4.2/2000.0;
    int order = 10;
    Teuchos::RCP<NRL_RandomFieldModel> interface = Teuchos::rcp(new NRL_RandomFieldModel(Comm,*paramList));
    
    double beta3 = -0.5; double beta4 = 2.495; double beta5 = 0.694; double plyagl = 0.26;
    interface->set_exponents(beta3,beta4,beta4);
    interface->set_plyagl(plyagl);
    interface->random_field_generator(seeds,mean,delta,corrlength,order);
    
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*interface,*paramList));
    Newton->Initialization();
    Newton->setParameters(*paramList);
    
    int error = Newton->Solve_with_Aztec(true);
    std::string filedisp   = "/Users/Brian/Documents/Thesis/Trilinos_results/nrl/randomfield/disp_realization3.mtx";
    std::string filestress = "/Users/Brian/Documents/Thesis/Trilinos_results/nrl/randomfield/stress_realization0.mtx";
    Newton->print_newton_solution(filedisp);
    //interface->compute_mean_cauchy_stress(*Newton->x,filename2);
        
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
