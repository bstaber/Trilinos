#ifndef LINEARFINITEELEMENTPROBLEM_HPP
#define LINEARFINITEELEMENTPROBLEM_HPP
#include "baseClassFEM.hpp"
class LinearFiniteElementProblem : public baseClassFEM
{
public:
    LinearFiniteElementProblem(){
    };
    ~LinearFiniteElementProblem(){
    };
    virtual void setup_dirichlet_conditions(){
    };
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
    };
};
#endif
