#ifndef NEARLYINCOMPRESSIBLEHYPERELASTICITY_HPP
#define NEARLYINCOMPRESSIBLEHYPERELASTICITY_HPP

#include "nearlyIncompressibleHyperelasticity.hpp"

class nearlyIncompressibleHyperelasticity : public hyperelasticity
{
public:
    nearlyIncompressibleHyperelasticity();
    ~nearlyIncompressibleHyperelasticity();
    
    void assemblePureDirichlet_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void assembleMixedDirichletDeformationDependentNeumann_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    
    void stiffnessRhsMaterialContribution(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    void stiffnessRhsPressureContribution(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
    
    virtual void get_constitutive_tensors_static_condensation(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseVector & inverse_cauchy, Epetra_SerialDenseVector & piola_isc, Epetra_SerialDenseVector & piola_vol, Epetra_SerialDenseMatrix & tangent_piola_isc, Epetra_SerialDenseMatrix & tangent_piola_vol) = 0;
    virtual void get_internal_pressure(double & theta, double & pressure, double & dpressure) = 0;
};
#endif
