/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef DAMAGEFIELD_HPP
#define DAMAGEFIELD_HPP

#include "linearFiniteElementProblem.hpp"

class damageField : public linearFiniteElementProblem
{
public:
    double gc;
    double lc;

    damageField(Epetra_Comm & comm, mesh & mesh, double & gc_, double & lc_);
    ~damageField();

    void solve(Teuchos::ParameterList & Parameters,
               Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
               Epetra_Vector & damageHistory, Epetra_Map & GaussMap);

    void assemble(Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
                  Epetra_Vector & damageHistory, Epetra_Map & GaussMap);

    void create_FECrsGraph();
    int print_solution(Epetra_Vector & lhs, std::string fileName);
    void setup_dirichlet_conditions();
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement);
};

#endif
