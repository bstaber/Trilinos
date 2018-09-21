#ifndef DAMAGEFIELD_HPP
#define DAMAGEFIELD_HPP

#include "linearFiniteElementProblem.hpp"

class damageField : public linearFiniteElementProblem
{
public:
    double gc;
    double lc;

    Epetra_Vector      * damageSolution;
    Epetra_FECrsMatrix * matrix;
    Epetra_FEVector    * rhs;

    damageField(Epetra_Comm & comm, Teuchos::ParameterList & Parameters, double & gc_, double & lc_);
    ~damageField();

    void solve(Teuchos::ParameterList & Parameters, Epetra_Vector & Hn);
    void assemble(Epetra_Vector & Hn);

    void create_FECrsGraph();
    void setup_dirichlet_conditions();
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement);
};

#endif
