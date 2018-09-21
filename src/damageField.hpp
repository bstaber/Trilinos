#ifndef DAMAGEFIELD_HPP
#define DAMAGEFIELD_HPP

include "linearFiniteElementProblem.hpp"

class damageField : public linearFiniteElementProblem
{
public:
    mesh * Mesh;
    double gc;
    double lc;

    Epetra_Vector damageSolution;

    damageField(mesh & Mesh, double & gc_, double & lc_);
    ~damageField();

    void solve();
    void assemble(Epetra_Vector & Hn, Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs);

    void create_FECrsGraph();
    void setup_dirichlet_conditions();
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement);

};

#endif
