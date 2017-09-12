#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Compressible_Mooney_Transverse_Isotropic.hpp"
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
    
    Epetra_SerialDenseVector parameters(8);
    Teuchos::RCP<TIMooney> interface = Teuchos::rcp(new TIMooney(Comm,*paramList));
    parameters(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
    parameters(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
    parameters(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
    parameters(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
    parameters(4) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
    parameters(5) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta3");
    parameters(6) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
    parameters(7) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
    double plyagl = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"angle");
    interface->set_parameters(parameters);
    interface->set_plyagl(plyagl);
    
    if (Comm.MyPID()==0){
    for (unsigned int i=0; i<parameters.Length(); i++){
        std::cout << parameters(i) << "\n";
    }
    }
    
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*interface,*paramList));
    
    std::vector<double> bcdisp(10);
    bcdisp[0] = 0.00033234/1000.0;
    bcdisp[1] = 0.018369/1000.0;
    bcdisp[2] = 0.038198/1000.0;
    bcdisp[3] = 0.060977/1000.0;
    bcdisp[4] = 0.073356/1000.0;
    bcdisp[5] = 0.092648/1000.0;
    bcdisp[6] = 0.11062/1000.0;
    bcdisp[7] = 0.12838/1000.0;
    bcdisp[8] = 0.14934/1000.0;
    bcdisp[9] = 0.0001571809118641;
    
    Newton->Initialization();
    for (unsigned int i=0; i<bcdisp.size(); ++i){
        Newton->setParameters(*paramList);
        Newton->bc_disp = bcdisp[i];
        int error = Newton->Solve_with_Aztec(true);
        //std::string name = "/Users/Brian/Documents/Thesis/Trilinos_results/nrl/deterministic/Newton_solution.mtx";
        //Newton->print_newton_solution(name);
    }
    
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
