#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "forwardDeterministicCostFunction.hpp"

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

    /*x(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
    x(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
    x(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
    x(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
    x(4) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
    x(5) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
    x(6) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
    for (unsigned int i=0; i<5; i++){
        x(i) = 1.0e3*x(i);
    }*/
    Epetra_SerialDenseVector x(7);
    /*x(0) = 3.208084e+03;
    x(1) = 4.516888e+01;
    x(2) = 4.024262e+01;
    x(3) = 1.587065e+03;
    x(4) = 2.603703e+01;
    x(5) = 1.335553e+01;
    x(6) = 3.304965e-02;*/

    /*x(0) = 2.739284e+03;
    x(1) = 2.613957e+01;
    x(2) = 4.524243e+01;
    x(3) = 1.785786e+03;
    x(4) = 3.367407e+01;
    x(5) = 2.527997e+01;
    x(6) = 4.677232e-02;*/

    x(0) = 1.734643e+03;
    x(1) = 6.033122e+01;
    x(2) = 3.537519e+01;
    x(3) = 1.265154e+03;
    x(4) = 4.468161e+01;
    x(5) = 2.426736e+01;
    x(6) = 3.924508e-02;

    {
        Teuchos::RCP<forwardDeterministicCostFunction> obj =
        Teuchos::rcp(new forwardDeterministicCostFunction(Comm,*paramList));
        Epetra_SerialDenseVector value(10);
        int id = 0;
        value = obj->value_id(x,id);
        if (Comm.MyPID()==0){
          std::cout << "************ ID = " << id << " ************\n";
          std::cout << value << "\n";
          std::cout << "\n";
        }
    }

    {
        Teuchos::RCP<forwardDeterministicCostFunction> obj =
        Teuchos::rcp(new forwardDeterministicCostFunction(Comm,*paramList));
        Epetra_SerialDenseVector value(10);
        int id = 1;
        value = obj->value_id(x,id);
        if (Comm.MyPID()==0){
          std::cout << "************ ID = " << id << " ************\n";
          std::cout << value << "\n";
          std::cout << "\n";
        }
    }

    {
        Teuchos::RCP<forwardDeterministicCostFunction> obj =
        Teuchos::rcp(new forwardDeterministicCostFunction(Comm,*paramList));
        Epetra_SerialDenseVector value(10);
        int id = 4;
        value = obj->value_id(x,id);
        if (Comm.MyPID()==0){
          std::cout << "************ ID = " << id << " ************\n";
          std::cout << value << "\n";
          std::cout << "\n";
        }
    }

    {
        Teuchos::RCP<forwardDeterministicCostFunction> obj =
        Teuchos::rcp(new forwardDeterministicCostFunction(Comm,*paramList));
        Epetra_SerialDenseVector value(10);
        int id = 6;
        value = obj->value_id(x,id);
        if (Comm.MyPID()==0){
          std::cout << "************ ID = " << id << " ************\n";
          std::cout << value << "\n";
          std::cout << "\n";
        }
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}