#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "forwardCostfunction.hpp"

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
    
    Teuchos::RCP<forwardCostfunction> costFunction = Teuchos::rcp(new forwardCostfunction(Comm,*paramList));
    Teuchos::RCP<readnrldata> data = Teuchos::rcp(new readnrldata(false));
    data->import_boundaryconditions();
    
    Epetra_IntSerialDenseVector seeds(5);
    int j = 0;
    for (int k=0; k<5; ++k){
        seeds(k) = 5*j+k;
    }
    
    int id = 0;
    Epetra_SerialDenseVector parameters(6), exponents(2), hyperParameters(7);
    hyperParameters(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"delta1");
    hyperParameters(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"delta2");
    hyperParameters(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"delta3");
    hyperParameters(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"delta4");
    hyperParameters(4) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"lx");
    hyperParameters(5) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"ly");
    parameters(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
    parameters(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
    parameters(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
    parameters(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
    parameters(4) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
    exponents(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
    exponents(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
    for (unsigned int i=0; i<5; i++){
        parameters(i) = 1.0e3*parameters(i);
    }

    //double value = costFunction->value(parameters,exponents,hyperParameters,id,seeds);
    
    /*double xi = 0.0;
    int j = 0;
    int * seed = new int [5];
    seed[0] = 5*j; seed[1] = 5*j+1; seed[2] = 5*j+2; seed[3] = 5*j+3; seed[4] = 5*j+4;
    interface->RandomFieldGenerator(seed);
    
    Newton->Initialization();
    for (unsigned int i=0; i<data->boundaryconditions.Length(); ++i){
        Newton->setParameters(*paramList);
        Newton->bc_disp = data->boundaryconditions(i);
        int error = Newton->Solve_with_Aztec(true);
        std::string path1 = "/home/s/staber/Trilinos_results/nrl/forward_randomfield/displacement_" + std::to_string(i) + ".mtx";
        std::string path2 = "/home/s/staber/Trilinos_results/nrl/forward_randomfield/greenlag_" + std::to_string(i) + ".mtx";
        Newton->print_newton_solution(path1);
        interface->compute_green_lagrange(*Newton->x,xi,xi,xi,path2);
    }*/
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
