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
    plasticitySmallStrains(Epetra_Comm & comm, Teuchos::ParameterList & parameterlist);
    ~plasticitySmallStrains();

    //void initialize(Teuchos::ParameterList & parameterlist);
    int incremental_bvp(bool print);

    void get_matrix_and_rhs(Epetra_Vector & eto, Epetra_Vector & ep_old, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    virtual void setup_dirichlet_conditions() = 0;
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement) = 0;

    void create_FECrsGraph();
    void constructGaussMap(Epetra_Map & Gauss_Map);

    Epetra_Map * GaussMap;

    int MyPID, iter_max, iter_min, nb_bis_max;
    double time, Delta, norm_inf_tol, norm_inf_max, eps, success, failure, umax, bc_disp, pressure_load, tol;

    unsigned int n_bc_dof;
    int * dof_on_boundary;

    Teuchos::ParameterList * Krylov;
};
#endif
