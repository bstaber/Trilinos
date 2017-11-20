#ifndef COMPRESSIBLEHYPERELASTICITY_HPP
#define COMPRESSIBLEHYPERELASTICITY_HPP

#include "hyperelasticity.hpp"

class compressibleHyperelasticity : public hyperelasticity
{
public:
    compressibleHyperelasticity();
    ~compressibleHyperelasticity();
    
    void assemblePureDirichlet_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void assemblePureDirichlet_inhomogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void assembleMixedDirichletNeumann_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void assembleMixedDirichletNeumann_inhomogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    
    void stiffnessRhs_homogeneousForcing(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void stiffnessRhs_inhomogeneousForcing(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void rhs_NeumannBoundaryCondition(Epetra_FEVector & F);
    
    virtual void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola) = 0;
    virtual Epetra_SerialDenseVector get_neumannBc(Epetra_SerialDenseMatrix & matrix_X, Epetra_SerialDenseMatrix & xg, unsigned int & gp) = 0;
    virtual Epetra_SerialDenseVector get_forcing(double & x1, double & x2, double & x3, unsigned int & e_lid, unsigned int & gp) = 0;
};
#endif
