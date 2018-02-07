#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Newton_Raphsonpp.hpp"
#include "RandomGeneratorForPCA_andLikelihood.hpp"
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
    if (Comm.MyPID()==0){
        paramList->print(std::cout,2,true,true);
    }

      Teuchos::RCP<RandomGeneratorForPCA_andLikelihood> RG =
      Teuchos::rcp(new RandomGeneratorForPCA_andLikelihood(Comm,*paramList));

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
      mean_parameters.Scale(1.0e3);
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
      plyagls(0) = 15.0*2.0*M_PI/360.0;
      plyagls(1) = 30.0*2.0*M_PI/360.0;
      plyagls(2) = 60.0*2.0*M_PI/360.0;
      plyagls(3) = 75.0*2.0*M_PI/360.0;

      int nmc = Teuchos::getParameter<int>(paramList->sublist("Shinozuka"),"nmc");
      Epetra_SerialDenseVector QoI(RG->nrldata->boundaryconditions.Length());
      Epetra_SerialDenseMatrix Z(RG->nrldata->boundaryconditions.Length(),4*nmc);

      int k = -1;
      for (unsigned int i=0; i<4; ++i){
        for (unsigned int j=0; j<nmc; ++j){
          k++;
          seeds(0) = 5*k+0;
          seeds(1) = 5*k+1;
          seeds(2) = 5*k+2;
          seeds(3) = 5*k+3;
          seeds(4) = 5*k+4;
          QoI = RG->rnd(seeds,
                        mean_parameters,
                        exponents,
                        correlation_lengths,
                        coeff_of_variation,
                        plyagls(i),
                        true,
                        false,
                        false);
          for (unsigned int l=0; l<QoI.Length(); ++l){
            //Z(l,j+i*nmc) = std::log(QoI(l));
          }
        }
      }

      if (Comm.MyPID()==0){
        std::ofstream output("/home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likelihood/output.txt");
        if (output.is_open()){
          output << Z;
          output.close();
        }
      }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
}
