/*
Brian Staber (brian.staber@gmail.com)
*/

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "tresca_plate.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

int main(int argc, char *argv[]){

    std::string xmlInFileName = "";

    Teuchos::CommandLineProcessor clp(false);
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
    Teuchos::RCP<tresca_plate> problem = Teuchos::rcp(new tresca_plate(Comm,*paramList));

    Epetra_SerialDenseVector DETO(6);
    //RIGHT, FIRST ITERATION
    DETO(0) = -0.000585561;
    DETO(1) = 0.0100571;
    DETO(2) = -0.000888424;
    DETO(3) = 0.000257445; //23
    DETO(4) = -2.42591e-06; //13
    DETO(5) = 0.00343659; //12
    Epetra_SerialDenseVector SIG(6);
    for (unsigned int k=0; k<6; ++k) SIG(k) = 0.0;
    double EPCUM = 0.0;

    //RIGHT, SECOND ITERATION
    //DETO(0) = 0.000407972; DETO(1) = 0.00424136; DETO(2) = -0.00135479; DETO(5) = -0.00243745; DETO(3) = 0.000392587; DETO(4) = -1.57293e-06;
    //SIG(0) = 1203.9; SIG(1) = 2144.82; SIG(2) = 1157.44; SIG(5) = 303.741; SIG(3) = 22.1086; SIG(4) = 4.80793;
    //EPCUM = 0.0;
    //EPCUM = 0.00345772; what it should be in Zset

    //SMOOTH, SECOND ITERATION
    /*DETO(0) = 0.000407972;
    DETO(1) = 0.00424136;
    DETO(2) = -0.00135479;
    DETO(5) = -0.00243745;
    DETO(3) = 0.000392587;
    DETO(4) = -1.57293e-06;
    SIG(0) = 1203.9; SIG(1) = 2144.82; SIG(2) = 1157.44; SIG(5) = 303.741; SIG(3) = 22.1086; SIG(4) = 4.80793;*/

    Epetra_SerialDenseMatrix TGM(6,6);

    unsigned int elid, gp;
    problem->constitutive_problem(elid,gp,DETO,SIG,EPCUM,TGM);

    Epetra_SerialDenseVector EELTEST(6);
    EELTEST.Multiply('N','N',1.0,problem->COMPLIANCE,SIG,0.0);

    std::cout << "DETO = " << DETO << std::endl;
    std::cout << "(LEFT,SMOOTH,RIGHT) = " << problem->LEFT << "\t" << problem->SMOOTH << "\t" << problem->RIGHT << std::endl;
    std::cout << "EEL_NEW = " << EELTEST << std::endl;
    std::cout << "SIG_NEW = " << SIG << std::endl;
    std::cout << "TGM_NEW = " << TGM << std::endl;
    std::cout << "EPCUM_NEW = " << EPCUM << std::endl;
    std::cout << "ELASTICITY = " << problem->ELASTICITY << std::endl;
    std::cout << "COMPLIANCE = " << problem->COMPLIANCE << std::endl;
    std::cout << "radius = " << problem->R0 + problem->H * EPCUM << std::endl;
    std::cout << problem->qs << "\t" << problem->qlsl << "\t" << problem->qlla << "\t" << problem->qrsr << "\t" << problem->qrra << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
