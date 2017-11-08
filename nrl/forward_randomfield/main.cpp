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
    
    int id = 0;
    Epetra_SerialDenseVector parameters(5), exponents(2), hyperParameters(6);
    parameters(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
    parameters(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
    parameters(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
    parameters(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
    parameters(4) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
    exponents(0)  = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
    exponents(1)  = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
    for (unsigned int i=0; i<5; i++){
        parameters(i) = 1.0e9*parameters(i);
    }
    
    if (Comm.MyPID()==0){
        std::cout << "\n";
        std::cout << "Value" << std::setw(10) << "Delta" << std::setw(10) << "lx" << std::setw(10) << "ly" << "\n";
    }
    
    double length = 50.0/1000.0;
    double width  = 25.0/1000.0;
    for (int I=1; I<=5; ++I){
        for (int J=1; J<=5; ++J){
            hyperParameters(0) = I/10.0;
            hyperParameters(1) = I/10.0;
            hyperParameters(2) = I/10.0;
            hyperParameters(3) = I/10.0;
            hyperParameters(4) = length*J*0.05;
            hyperParameters(5) = width*J*0.05;
            for (int j=0; j<10; ++j){
                for (int k=0; k<5; ++k){
                    seeds(k) = 5*j+k;
                }
                double value = costFunction->value(parameters,exponents,hyperParameters,id,seeds,false);
                std::string path1 = "/home/s/staber/Trilinos_results/nrl/forward_randomfield/u_delta" + std::to_string(I) + "_L" + std::to_string(J) + "_nmc" + std::to_string(j) + ".mtx";
                std::string path2 = "/home/s/staber/Trilinos_results/nrl/forward_randomfield/e_delta" + std::to_string(I) + "_L" + std::to_string(J) + "_nmc" + std::to_string(j) + ".mtx";
                costFunction->print_newton_solution(path1);
                costFunction->print_green_lagrange(path2);
                if (Comm.MyPID()==0){
                    std::cout << value << std::setw(10) << I/10.0 << std::setw(10) << length*J*0.05 << std::setw(10) << width*J*0.05 << "\n";
                }
            }
        }
    }
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
