#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "objectiveFunction.hpp"
#include <random>

int main(int argc, char *argv[]){
    
    Teuchos::RCP<objectiveFunction<double>> obj = Teuchos::rcp(new objectiveFunction<double>());
    
    std::string path1 = "/Users/Brian/Documents/Thesis/0-Trilinos/Trilinos/Fitting/data/specimen1_circ.txt";
    std::string path2 = "/Users/Brian/Documents/Thesis/0-Trilinos/Trilinos/Fitting/data/specimen1_long.txt";
    
    obj->set_experimental_data(path1,path2);
    
    Epetra_SerialDenseVector p(7);
    p(0) = 3.0759;
    p(1) = 0.21828;
    p(2) = 5.2622;
    p(3) = 37.259;
    p(4) = 3.5384;
    p(5) = 496.49;
    p(6) = 0.82323;
    
    double tau;
    Epetra_SerialDenseVector y(3);
    
    std::mt19937 G(time(NULL));
    std::uniform_real_distribution<double>  mu4(1.0,1000.0);
    std::uniform_real_distribution<double>  beta4(2.2,1000.0);
    
    for (unsigned int nmc=0; nmc<1e6; nmc++){
        
        p(3) = mu4(G);
        p(5) = beta4(G);
        
        y(1) = 1.0;
        y(2) = 1.0;
    
        std::cout << "\n";
        std::cout << std::setw(15) << "residual" << std::setw(10) << "y1" << std::setw(10) << "y2" << std::setw(10) << "y3" << std::setw(10) << "tau" << std::setw(10) << "mu1 " << std::setw(10) << "mu2" << std::setw(10) << "mu3" << std::setw(10) << "mu4" << std::setw(10) << "beta3" << std::setw(10) << "beta4" << std::setw(10) << "phi" << "\n";
        std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
        
        for (unsigned int i=0; i<obj->y1_exp_circ.Length(); ++i){
            y(0) = obj->y1_exp_circ(i);
            double fval = obj->newton(y,p);
            obj->compute_stress(y,p,tau);
            std::cout << std::setw(15) << fval << std::setw(10) << y(0) << std::setw(10) << y(1) << std::setw(10) << y(2) << std::setw(10) << tau;
            for (unsigned int k=0; k<p.Length(); k++){
                std::cout << std::setw(10) << p(k);
            }
            std::cout << "\n";
        }
        
    }
    
    return 0;
}
