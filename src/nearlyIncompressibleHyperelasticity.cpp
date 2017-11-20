#include "nearlyIncompressibleHyperelasticity.hpp"

nearlyIncompressibleHyperelasticity::nearlyIncompressibleHyperelasticity(){
}

nearlyIncompressibleHyperelasticity::~nearlyIncompressibleHyperelasticity(){
}

void nearlyIncompressibleHyperelasticity::assemblePureDirichlet_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
void nearlyIncompressibleHyperelasticity::assembleMixedDirichletDeformationDependentNeumann_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F);

void nearlyIncompressibleHyperelasticity::stiffnessRhsMaterialContribution(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
void nearlyIncompressibleHyperelasticity::stiffnessRhsPressureContribution(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F);
