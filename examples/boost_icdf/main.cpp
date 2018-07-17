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

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

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

    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> phi_(0.0,0.1);

    rng.seed(0);
    std::cout << phi_(rng) << "\n";
    std::cout << phi_(rng) << "\n";
    std::cout << phi_(rng) << "\n";
    rng.seed(0);
    std::cout << phi_(rng) << "\n";
    std::cout << phi_(rng) << "\n";
    std::cout << phi_(rng) << "\n";
    /*double mean = 0.0;
    double stan = 1.0;
    boost::random::normal_distribution<> w(mean,stan);
    boost::random::mt19937 rng;
    rng.seed(std::time(0));
    for (unsigned int i=0; i<10; ++i){
        double x = w(rng);
        double erfarg = (x-mean)/(stan*std::sqrt(2.0));
        double erfx = boost::math::erf<double>(erfarg);
        double y = (1.0/2.0)*(1.0 + erfx);
        double yinv = boost::math::gamma_p_inv<double,double>(1.0/(0.1*0.1),y);
        double z = yinv*10.0*0.1*0.1;
        double z2 = boost::math::ibeta_inv<double,double,double>(5.0,4.0,y);
        std::cout << x << std::setw(20) << y << std::setw(20) << z << std::setw(20) << z2 << "\n";
    }*/


#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
