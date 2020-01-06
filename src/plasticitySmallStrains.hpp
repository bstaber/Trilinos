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

    void initialize(Epetra_Comm & comm, Teuchos::ParameterList & parameterlist);

    int incremental_bvp(bool print);
    int sequence_bvp(bool print);

    void assemble_system(const Epetra_Vector & Du_, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void stiffness_rhs_homogeneousForcing(const Epetra_Vector & Du_, Epetra_FECrsMatrix & K, Epetra_FEVector & F);

    void elastic_predictor(Epetra_LinearProblem & problem_, AztecOO & solver_, Epetra_FECrsMatrix & K, Epetra_FEVector & rhs_,
                           Epetra_Vector & lhs_, double & displacement_);
    void assemble_system_LinearElasticity(Epetra_FECrsMatrix & K);
    void stiffness_homogeneousForcing_LinearElasticity(Epetra_FECrsMatrix & K);

    void compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B);

    void solve(Epetra_LinearProblem & problem_, AztecOO & solver_,
               Epetra_FECrsMatrix & A, Epetra_FEVector & b, Epetra_Vector & lhs_);
    void integrate_constitutive_problem(Epetra_Vector & Du_);

    virtual void constitutive_problem(const unsigned int & elid, const unsigned int & igp,
                                      const Epetra_SerialDenseVector & DETO, Epetra_SerialDenseVector & SIG,
                                      double & EPCUM, Epetra_SerialDenseMatrix & TGM) = 0;
    virtual void get_elasticity_tensor(const unsigned int & elid, const unsigned int & igp, Epetra_SerialDenseMatrix & TGM) = 0;

    virtual void setup_dirichlet_conditions() = 0;
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double bc) = 0;

    void create_FECrsGraph();
    void constructScalarGaussMap();

    int print_solution(Epetra_Vector & solution, std::string fileName);
    int print_at_gauss_points(std::string filename);
    int print_mean_at_center(std::string filename);

    Teuchos::ParameterList * Krylov;
    Teuchos::ParameterList * Newton;

    Epetra_Time * CpuTime;

    int MyPID, iter_max, iter_min, nb_bis_max;
    double time, Delta, norm_inf_tol, norm_inf_max, eps, success, failure, umax, bc_disp, pressure_load, tol;
    double Assemble_time, Aztec_time;

    unsigned int n_bc_dof;
    int * dof_on_boundary;

    Epetra_Map         * GaussMap;
    Epetra_Vector      * epcum, * epcum_converged;
    Epetra_MultiVector * sig, * sig_converged;
    Epetra_MultiVector * eto, * eto_converged;
    Epetra_MultiVector * tgm;

    Epetra_SerialDenseMatrix ELASTICITY;
};
#endif
