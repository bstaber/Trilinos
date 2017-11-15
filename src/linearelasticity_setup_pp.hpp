#ifndef LINEARELASTICITY_SETUP_PP_HPP
#define LINEARELASTICITY_SETUP_PP_HPP

#include "Linear_Finite_Element_Problem.hpp"

class LinearizedElasticity : public Linear_Finite_Element_Problem
{
public:
    LinearizedElasticity();
    ~LinearizedElasticity();
    
    void create_FECrsGraph();
    
    void assemblePureDirichlet_homogeneousForcing(Epetra_FECrsMatrix & K);
    void assembleMixedDirichletNeumann_homogeneousForcing(Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    
    void stiffness_pureDirichlet_homogeneousForcing(Epetra_FECrsMatrix & K);
    void stiffness_pureDirichlet_inhomogeneousForcing(Epetra_FECrsMatrix & K, Epetra_FECrsMatrix & F);
    void rhs_NeumannBoundaryCondition(Epetra_FEVector & F);
    
    void compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B);
    
    void compute_mean_cauchy_stress(Epetra_Vector & x, std::string & filename, bool printCauchy, bool printVM);
    void compute_deformation(Epetra_Vector & x, std::string & filename, bool printCauchy, bool printVM);
    
    virtual Epetra_SerialDenseVector get_neumannBc(unsigned int & e_lid, unsigned int & gp) = 0;
    virtual Epetra_SerialDenseVector get_forcing(unsigned int & e_lid, unsigned int & gp) = 0;
    virtual void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix) = 0;
    virtual void get_elasticity_tensor_for_recovery(unsigned int & e_lid, Epetra_SerialDenseMatrix & tangent_matrix) = 0;
    
    unsigned int n_bc_dof;
    int * dof_on_boundary;    
};
#endif
