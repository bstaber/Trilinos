#ifndef HYPERELASTICITY_HPP
#define HYPERELASTICITY_HPP

#include "Finite_Element_Problem.hpp"

class hyperelasticity : public Finite_Element_Problem
{
public:
    hyperelasticity();
    ~hyperelasticity();

    void create_FECrsGraph();
    void compute_B_matrices(Epetra_SerialDenseMatrix & F, Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B, Epetra_SerialDenseMatrix & BG);

    void compute_green_lagrange(Epetra_Vector & x, double & xi, double & eta, double & zeta, std::string & filename);
    void compute_center_cauchy_stress(Epetra_Vector & x, std::string & filename);
    void compute_gauss_vonmises(Epetra_Vector & x, std::string & filename);

    unsigned int n_bc_dof;
    int * dof_on_boundary;

    virtual void get_material_parameters(unsigned int & e_lid, unsigned int & gp) = 0;
    virtual void get_material_parameters_for_recover(unsigned int & e_lid) = 0;
    virtual void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress) = 0;
};
#endif
