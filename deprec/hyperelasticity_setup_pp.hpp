#ifndef HYPERELASTICITY_SETUP_PP_HPP
#define HYPERELASTICITY_SETUP_PP_HPP

#include "Finite_Element_Problem.hpp"
#include "laplacepp.hpp"

class hyperelasticity_setup : public Finite_Element_Problem
{
public:
    hyperelasticity_setup();
    ~hyperelasticity_setup();
    
    void create_FECrsGraph();
    
    void assemble_dirichlet_live_neumann_static_condensation(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void assemble_dirichlet_static_condensation(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    
    void material_stiffness_and_rhs_static_condensation(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void force_stiffness_rhs_live_pressure(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
        
    void compute_B_matrices(Epetra_SerialDenseMatrix & F, Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B, Epetra_SerialDenseMatrix & BG);
    
    void compute_green_lagrange(Epetra_Vector & x, double & xi, double & eta, double & zeta, std::string & filename);
    void compute_mean_cauchy_stress(Epetra_Vector & x, std::string & filename);
    void compute_gauss_vonmises(Epetra_Vector & x, std::string & filename);
    
    virtual void get_material_parameters(unsigned int & e_lid, unsigned int & gp) = 0;
    virtual void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola) = 0;
    virtual void get_constitutive_tensors_static_condensation(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseVector & inverse_cauchy, Epetra_SerialDenseVector & piola_isc, Epetra_SerialDenseVector & piola_vol, Epetra_SerialDenseMatrix & tangent_piola_isc, Epetra_SerialDenseMatrix & tangent_piola_vol) = 0;
    virtual void get_internal_pressure(double & theta, double & pressure, double & dpressure) = 0;
    virtual void get_material_parameters_for_recover(unsigned int & e_lid, double & xi, double & eta, double & zeta) = 0;
    virtual void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress) = 0;
    
    unsigned int n_bc_dof;
    int * dof_on_boundary;    
};
#endif
