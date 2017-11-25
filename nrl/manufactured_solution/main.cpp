#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "manufacturedSolution.hpp"

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
    if(xmlInFileName.length()){
        Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(*paramList));
    }
    
    if (Comm.MyPID()==0){
        paramList->print(std::cout,2,true,true);
    }
    
    Epetra_SerialDenseVector parameters(7);
    parameters(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
    parameters(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
    parameters(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
    parameters(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
    parameters(4) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
    parameters(5) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
    parameters(6) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
    for (unsigned int i=0; i<5; ++i){
        parameters(i) = 1.0e3*parameters(i);
    }
    double plyagl = 2.0*M_PI*30.0/360.0;
    
    unsigned int n = 5;
    Epetra_SerialDenseVector errorL2(n);
    for (unsigned int i=0; i<n; ++i){
        //std::string mesh_file  = "/Users/brian/Documents/GitHub/Trilinos/cee530/mesh/manufactured" + std::to_string(i) + ".msh";
        std::string mesh_file  = "/home/s/staber/Trilinos/cee530/mesh/manufactured" + std::to_string(i) + ".msh";
    
        Teuchos::RCP<manufacturedSolution> manufactured = Teuchos::rcp(new manufacturedSolution(Comm,*paramList,mesh_file));
        Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*manufactured,*paramList));
        manufactured->set_parameters(parameters,plyagl);
    
        Newton->Initialization();
        Newton->setParameters(*paramList);
        int error = Newton->Solve_with_Aztec(true);
        //std::string path1 = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/manufactured/manufactured" + std::to_string(i) + ".mtx";
        std::string path1 = "/home/s/staber/Trilinos_results/nrl/manufactured/manufactured" + std::to_string(i) + ".mtx";
        Newton->print_newton_solution(path1);
        errorL2(i) = manufactured->errorL2(*Newton->x);
    }
    Comm.Barrier();
    if (Comm.MyPID()==0){
        std::cout << errorL2;
    }
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}

//0.005 0.0105263157894737 0.0160526315789474 0.0215789473684211 0.0271052631578947 0.0326315789473684 0.0381578947368421 0.0436842105263158 0.0492105263157895 0.0547368421052632 0.0602631578947368 0.06578947368421049 0.07131578947368419 0.07684210526315791 0.0823684210526316 0.0878947368421053 0.0934210526315789 0.0989473684210526 0.1044736842105263 0.11 
