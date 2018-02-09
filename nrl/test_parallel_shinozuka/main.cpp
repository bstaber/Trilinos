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
#include "Compressible_Mooney_Transverse_Isotropic_Random_Field.hpp"

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

    Teuchos::RCP<TIMooney_RandomField> interface =
    Teuchos::rcp(new TIMooney_RandomField(Comm,*paramList));

    Epetra_IntSerialDenseVector seeds(5);
    Epetra_SerialDenseVector    mean_parameters(5);
    Epetra_SerialDenseVector    exponents(2);
    Epetra_SerialDenseVector    omega(6);
    Epetra_SerialDenseVector    plyagls(4);

    //mean values of the random parameters G_1(x),...,G_5(x)
    mean_parameters(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
    mean_parameters(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
    mean_parameters(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
    mean_parameters(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
    mean_parameters(4) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
    mean_parameters.Scale(1.0e3);
    //deterministic exponents beta_4 and beta_5
    exponents(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
    exponents(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
    //correlation lengths of the Gaussian random field
    omega(4) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"lx");
    omega(5) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"ly");
    //coefficients of variation of the random parameters G_1(x),...,G_4(x)
    omega(0) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"delta1");
    omega(1) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"delta2");
    omega(2) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"delta3");
    omega(3) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"delta4");
    //ply angle
    plyagls(0) = 15.0*2.0*M_PI/360.0;
    plyagls(1) = 30.0*2.0*M_PI/360.0;
    plyagls(2) = 60.0*2.0*M_PI/360.0;
    plyagls(3) = 75.0*2.0*M_PI/360.0;

    for (unsigned int i=0; i<5; ++i){ seeds(i) = i;};
    interface->setParameters(mean_parameters,exponents,omega);
    interface->set_plyagl(plyagls(0));
    interface->RandomFieldGenerator(seeds);

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
