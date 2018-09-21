#ifndef PHASEFIELDLINEARIZEDELASTICITY_HPP
#define PHASEFIELDLINEARIZEDELASTICITY_HPP

#include "linearizedElasticity.hpp"

class phaseFieldLinearizedElasticity : public linearizedElasticity{

  Epetra_SerialDenseMatrix elasticity_tensor;
  Teuchos::RCP<damageField> damageInterface;
  double gc;
  double lc;

  Epetra_Vector * damageHistory;

  phaseFieldLinearizedElasticity(mesh & mesh, double & gc_, double & lc_);
  ~phaseFieldLinearizedElasticity();

  void computeDisplacement();
  void computeDamageHistory();

  void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix);
  void get_elasticity_tensor_for_recovery(unsigned int & e_lid, Epetra_SerialDenseMatrix & tangent_matrix);

};
#endif
