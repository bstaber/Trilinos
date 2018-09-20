#ifndef PHASEFIELDLINEARIZEDELASTICITY_HPP
#define PHASEFIELDLINEARIZEDELASTICITY_HPP

#include "linearizedElasticity.hpp"

class phaseFieldLinearizedElasticity : public linearizedElasticity{

  Teuchos::RCP<damageField> damage;
  double gc;
  double lc;

  phaseFieldLinearizedElasticity(mesh & mesh, double & gc_, double & lc_);
  ~phaseFieldLinearizedElasticity();

  void solve();

  Epetra_SerialDenseVector get_neumannBc(Epetra_SerialDenseMatrix & matrix_X, Epetra_SerialDenseMatrix & xg, unsigned int & gp);
  Epetra_SerialDenseVector get_forcing(double & x1, double & x2, double & x3, unsigned int & e_lid, unsigned int & gp);
  void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix);
  void get_elasticity_tensor_for_recovery(unsigned int & e_lid, Epetra_SerialDenseMatrix & tangent_matrix);
  void setup_dirichlet_conditions();
  void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement);

};
#endif
