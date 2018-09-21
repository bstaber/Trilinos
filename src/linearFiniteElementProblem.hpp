/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef LINEARFINITEELEMENTPROBLEM_HPP
#define LINEARFINITEELEMENTPROBLEM_HPP
#include "baseClassFEM.hpp"
class linearFiniteElementProblem : public baseClassFEM
{
public:
    linearFiniteElementProblem(){
    };
    ~linearFiniteElementProblem(){
    };
    virtual void setup_dirichlet_conditions(){
    };
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
    };
};
#endif
