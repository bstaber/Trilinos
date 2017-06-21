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

#include "meshpp.hpp"

int main(int argc, char *argv[]){
    
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    
    double mean = 0.0;
    double stan = 1.0;
    int nmc = 1000;
    int n = 3;
    
    boost::random::normal_distribution<double> x(mean,stan);
    boost::random::mt19937 rng;
    rng.seed(std::time(0));
    
    Epetra_SerialDenseVector X(n);
    Epetra_SerialDenseVector Y(n);
    Epetra_SerialDenseVector c(n);
    Epetra_SerialDenseMatrix Z(nmc,n);
    c(0) = 0.0; c(1) = 0.0; c(2) = 0.0;
    double radius = 1.0;
    
    for (unsigned int i=0; i<nmc; ++i){
        for (unsigned int j=0; j<n; ++j){
            X(j) = x(rng);
        }
        double s = X.Norm2();
        double gammainc = boost::math::gamma_p<double,double>(s*s/2.0,double(n)/2.0);
        for (unsigned int j=0; j<n; ++j){
            Y(j) = X(j)*(radius*std::pow(gammainc,1.0/double(n))/s);
            Z(i,j) = Y(j) + c(j);
            std::cout << std::setw(20) << Z(i,j);
        }
        std::cout << "\n";
    }
    
 
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
