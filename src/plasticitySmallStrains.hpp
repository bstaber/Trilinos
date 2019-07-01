/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef PLASTICITYSMALLSTRAINS_HPP
#define PLASTICITYSMALLSTRAINS_HPP

#include "nonLinearFiniteElementProblem.hpp"
#include "newtonRaphson.hpp"

class plasticitySmallStrains : public baseClassFEM
{
public:
    plasticitySmallStrains();
    ~plasticitySmallStrains();

    void create_FECrsGraph();

    void incremental_problem();
    void global_newton();

    void get_matrix_and_rhs(Epetra_Vector & eto, Epetra_Vector & ep_old, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    virtual void setup_dirichlet_conditions() = 0;
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement) = 0;

    unsigned int n_bc_dof;
    int * dof_on_boundary;
};
#endif
