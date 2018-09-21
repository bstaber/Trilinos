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

#include "nrl_PCA_Likelihood.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>

int main(int argc, char *argv[]){

    std::string xmlInFileName = "";

    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setDocString("Compilation: make -f Makefile.os\n"
                     "Run: mpirun -np 28 ./trilinos_mpi --xml-in-file='nrl.aztec.linux.xml'");

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

    /*if (Comm.MyPID()==0){
        std::cout << "****************\n";
        std::cout << "Hyperparameters:\n";
        std::cout << paramList->sublist("TIMooney");
        std::cout << paramList->sublist("Shinozuka");
        std::cout << "****************\n";
    }*/

      Teuchos::RCP<nrl_PCA_Likelihood> RG =
      Teuchos::rcp(new nrl_PCA_Likelihood(Comm,*paramList));

      Epetra_IntSerialDenseVector seeds(5);
      Epetra_SerialDenseVector    mean_parameters(5);
      Epetra_SerialDenseVector    exponents(2);
      Epetra_SerialDenseVector    correlation_lengths(2);
      Epetra_SerialDenseVector    coeff_of_variation(4);
      Epetra_SerialDenseVector    plyagls(4);

      //mean values of the random parameters G_1(x),...,G_5(x)
      mean_parameters(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
      mean_parameters(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
      mean_parameters(2) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
      mean_parameters(3) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
      mean_parameters(4) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
      //mean_parameters.Scale(1.0e3);
      //deterministic exponents beta_4 and beta_5
      exponents(0) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
      exponents(1) = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
      //correlation lengths of the Gaussian random field
      correlation_lengths(0) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"lx");
      correlation_lengths(1) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"ly");
      //coefficients of variation of the random parameters G_1(x),...,G_4(x)
      coeff_of_variation(0) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"delta1");
      coeff_of_variation(1) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"delta2");
      coeff_of_variation(2) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"delta3");
      coeff_of_variation(3) = Teuchos::getParameter<double>(paramList->sublist("Shinozuka"),"delta4");
      //ply angle
      plyagls(0) = 15.0;
      plyagls(1) = 30.0;
      plyagls(2) = 60.0;
      plyagls(3) = 75.0;

      int nmc = Teuchos::getParameter<int>(paramList->sublist("Shinozuka"),"nmc");

      for (unsigned int i=0; i<4; ++i){
        for (unsigned int j=0; j<nmc; ++j){
          seeds(0) = 5*(j+i*nmc)+0;
          seeds(1) = 5*(j+i*nmc)+1;
          seeds(2) = 5*(j+i*nmc)+2;
          seeds(3) = 5*(j+i*nmc)+3;
          seeds(4) = 5*(j+i*nmc)+4;
          int flag = RG->rnd(j                  ,
                             seeds              ,
                             mean_parameters    ,
                             exponents          ,
                             correlation_lengths,
                             coeff_of_variation ,
                             plyagls(i)         ,
                             false              ,
                             false              ,
                             false              ,
                             true               );
        }
      }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
}
