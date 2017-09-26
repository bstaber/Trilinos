#ifndef LINEARELASTICITY_SETUP_PP_HPP
#define LINEARELASTICITY_SETUP_PP_HPP

#include "Linear_Finite_Element_Problem.hpp"

class LinearizedElasticity : public Linear_Finite_Element_Problem
{
public:
    LinearizedElasticity();
    ~LinearizedElasticity();
    
    void create_FECrsGraph();
    
    void assemble_dirichlet(Epetra_FECrsMatrix & K);
    void assemble_dirichlet_dead_neumann(Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    
    void material_stiffness_and_rhs_dirichlet(Epetra_FECrsMatrix & K);
    void force_dead_pressure(Epetra_FEVector & F);
    
    void compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B);
    
    void compute_mean_cauchy_stress(Epetra_Vector & x, std::string & filename, bool printCauchy, bool printVM);
    void compute_deformation(Epetra_Vector & x, std::string & filename, bool printCauchy, bool printVM);
    
    virtual void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix) = 0;
    virtual void get_elasticity_tensor_for_recovery(unsigned int & e_lid, Epetra_SerialDenseMatrix & tangent_matrix) = 0;
    
    unsigned int n_bc_dof;
    int * dof_on_boundary;    
};
#endif
