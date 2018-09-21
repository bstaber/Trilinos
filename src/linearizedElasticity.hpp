/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef LINEARIZEDELASTICITY_HPP
#define LINEARIZEDELASTICITY_HPP

#include "linearFiniteElementProblem.hpp"

class linearizedElasticity : public linearFiniteElementProblem
{
public:
    linearizedElasticity();
    ~linearizedElasticity();

    void create_FECrsGraph();

    void aztecSolver(Epetra_FECrsMatrix & A, Epetra_FEVector & b, Epetra_Vector & u, Teuchos::ParameterList & paramList);

    void assemblePureDirichlet_homogeneousForcing(Epetra_FECrsMatrix & K);
    void assembleMixedDirichletNeumann_homogeneousForcing(Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void assembleMixedDirichletNeumann_inhomogeneousForcing(Epetra_FECrsMatrix & K, Epetra_FEVector & F);

    void stiffness_homogeneousForcing(Epetra_FECrsMatrix & K);
    void stiffness_inhomogeneousForcing(Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void rhs_NeumannBoundaryCondition(Epetra_FEVector & F);

    void compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B);

    void compute_center_cauchy_stress(Epetra_Vector & x, std::string & filename, bool printCauchy, bool printVM);
    void compute_deformation(Epetra_Vector & x, std::string & filename, bool printCauchy, bool printVM);

    int print_solution(Epetra_Vector & solution, std::string fileName);

    virtual Epetra_SerialDenseVector get_neumannBc(Epetra_SerialDenseMatrix & matrix_X, Epetra_SerialDenseMatrix & xg, unsigned int & gp) = 0;
    virtual Epetra_SerialDenseVector get_forcing(double & x1, double & x2, double & x3, unsigned int & e_lid, unsigned int & gp) = 0;
    virtual void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix) = 0;
    virtual void get_elasticity_tensor_for_recovery(unsigned int & e_lid, Epetra_SerialDenseMatrix & tangent_matrix) = 0;

    unsigned int n_bc_dof;
    int * dof_on_boundary;
};
#endif
