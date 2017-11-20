#ifndef COMPRESSIBLEHYPERELASTICITY_HPP
#define COMPRESSIBLEHYPERELASTICITY_HPP

#include "compressibleHyperelasticity.hpp"

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
};
#endif
