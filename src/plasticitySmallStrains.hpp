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

    void assemblePureDirichlet_homogeneousForcing(const Epetra_Vector & du, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void assemblePureDirichlet_homogeneousForcing_LinearElasticity(Epetra_FECrsMatrix & K);
    void stiffness_rhs_homogeneousForcing(const Epetra_Vector & du, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void stiffness_homogeneousForcing_LinearElasticity(Epetra_FECrsMatrix & K);
    void compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B);
    virtual void constitutive_problem(const unsigned int & elid, const unsigned int & igp,
                              Epetra_SerialDenseVector & deto, Epetra_SerialDenseVector & sig,
                              Epetra_SerialDenseMatrix & tgm) = 0;
    virtual void get_elasticity_tensor(unsigned int & elid, unsigned int & igp, Epetra_SerialDenseMatrix & tgm) = 0;

    virtual void setup_dirichlet_conditions() = 0;
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement) = 0;

    void create_FECrsGraph();
    void constructGaussMap(Epetra_Map & Gauss_Map);

    Teuchos::ParameterList * Krylov;

    int MyPID, iter_max, iter_min, nb_bis_max;
    double time, Delta, norm_inf_tol, norm_inf_max, eps, success, failure, umax, bc_disp, pressure_load, tol;

    unsigned int n_bc_dof;
    int * dof_on_boundary;

    Epetra_Map * GaussMap;
    Epetra_Vector * eto;
    Epetra_Vector * eel;
    Epetra_Vector * epi;

    Epetra_SerialDenseMatrix ELASTICITY;
};
#endif
