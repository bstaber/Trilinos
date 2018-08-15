#ifndef NONLINEARFINITEELEMENTPROBLEM_HPP
#define NONLINEARFINITEELEMENTPROBLEM_HPP

#include "meshpp.hpp"
#include "baseClassFEM.hpp"

class nonLinearFiniteElementProblem : public baseClassFEM
{
    
public:
    nonLinearFiniteElementProblem(){
    };
    ~nonLinearFiniteElementProblem(){
    };
    
    double pressure_load;
    
    virtual void get_matrix_and_rhs(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F) = 0;
    virtual void setup_dirichlet_conditions() = 0;
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement) = 0;
    
    void display_amesos_solvers(){
        Amesos Factory;
        std::cout << "Available Amesos solvers: \n";
        bool isAvailable;
        isAvailable = Factory.Query("Lapack"); if(isAvailable){ std::cout << "Lapack: yes.\n"; } else{ std::cout << "Lapack: no.\n"; }
        isAvailable = Factory.Query("Klu"); if(isAvailable){ std::cout << "Klu: yes.\n"; } else{ std::cout << "Klu: no.\n"; }
        isAvailable = Factory.Query("Umfpack"); if(isAvailable){ std::cout << "Umfpack: yes.\n"; } else{ std::cout << "Umfpack: no.\n"; }
        isAvailable = Factory.Query("Pardiso"); if(isAvailable){ std::cout << "Pardiso: yes.\n"; } else{ std::cout << "Pardiso: no.\n"; }
        isAvailable = Factory.Query("Taucs"); if(isAvailable){ std::cout << "Taucs: yes.\n"; } else{ std::cout << "Taucs: no.\n"; }
        isAvailable = Factory.Query("Superlu"); if(isAvailable){ std::cout << "Superlu: yes.\n"; } else{ std::cout << "Superlu: no.\n"; }
        isAvailable = Factory.Query("Superludist"); if(isAvailable){ std::cout << "Superludist: yes.\n"; } else{ std::cout << "Superludist: no.\n"; }
        isAvailable = Factory.Query("Mumps"); if(isAvailable){ std::cout << "Mumps: yes.\n"; } else{ std::cout << "Mumps: no.\n"; }
        isAvailable = Factory.Query("Dscpack"); if(isAvailable){ std::cout << "Dscpack: yes.\n"; } else{ std::cout << "Dscpack: no.\n"; }
    }
    
};

#endif
