/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef PHASEFIELDPROBLEM_HPP
#define PHASEFIELDPROBLEM_HPP

#include "phaseFieldLinearizedElasticity.hpp"

class phaseFieldProblem : public phaseFieldLinearizedElasticity{

public:
  phaseFieldProblem(Epetra_Comm & comm, Teuchos::ParameterList & Parameters);
  ~phaseFieldProblem();

  Epetra_SerialDenseVector get_neumannBc(Epetra_SerialDenseMatrix & matrix_X, Epetra_SerialDenseMatrix & xg, unsigned int & gp);
  Epetra_SerialDenseVector get_forcing(double & x1, double & x2, double & x3, unsigned int & e_lid, unsigned int & gp);
  void setup_dirichlet_conditions();
  void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement);

};
#endif
