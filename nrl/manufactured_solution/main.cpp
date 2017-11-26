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
    int            meshIndex;
    
    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setOption("mesh-index",&meshIndex,"The index of the manufactured mesh.");
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
    
    double errorL2;
    std::string path     = Teuchos::getParameter<std::string>(paramList->sublist("Mesh"),"path");
    std::string meshname = "manufactured" + std::to_string(meshIndex) + ".msh";
    std::string fullpath = path + meshname;
    Teuchos::RCP<manufacturedSolution> manufactured = Teuchos::rcp(new manufacturedSolution(Comm,*paramList,fullpath));
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*manufactured,*paramList));
    manufactured->set_parameters(parameters,plyagl);
    Comm.Barrier();
        
    Newton->Initialization();
    Newton->setParameters(*paramList);
    int error = Newton->Solve_with_Aztec(true);
    
    std::string outpath     = Teuchos::getParameter<std::string>(paramList->sublist("Mesh"),"outpath");
    std::string fulloutpath = outpath + std::to_string(meshIndex) + ".mtx";
    Newton->print_newton_solution(fulloutpath);
    
    double xi = 0.0;
    fulloutpath = outpath + "_gl_" + std::to_string(meshIndex) + ".mtx";
    manufactured->compute_green_lagrange(*Newton->x,xi,xi,xi,fulloutpath);
    errorL2 = manufactured->errorL2(*Newton->x);
    
    if (Comm.MyPID()==0){
        std::ofstream output;
        output.open("output.txt",std::ios_base::app);
        output << meshname << std::setw(5) << Comm.NumProc() << std::setw(15) << errorL2 << "\n";
        output.close();
    }
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
