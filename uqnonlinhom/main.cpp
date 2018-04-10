#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "StochasticHomogenization.hpp"
#include "Newton_Raphsonpp.hpp"
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

      Teuchos::RCP<StochasticHomogenization> interface =
      Teuchos::rcp(new StochasticHomogenization(Comm,*paramList));

      Teuchos::RCP<Newton_Raphson> newton =
      Teuchos::rcp(new Newton_Raphson(*interface,*paramList));

      Epetra_SerialDenseVector x(6);
      x(0) = 1.0;
      x(1) = 1.0;
      x(2) = 1.0;
      x(3) = 1.0;
      x(4) = 1.0;
      x(5) = 1.0;
      interface->set_parameters(x);


#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
}
