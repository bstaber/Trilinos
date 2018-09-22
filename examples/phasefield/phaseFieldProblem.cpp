/*
Brian Staber (brian.staber@gmail.com)
*/


#include "phaseFieldProblem.hpp"

phaseFieldProblem::phaseFieldProblem(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
  initialize(comm, Parameters);
}

phaseFieldProblem::~phaseFieldProblem(){
}

Epetra_SerialDenseVector phaseFieldProblem::get_neumannBc(Epetra_SerialDenseMatrix & matrix_X, Epetra_SerialDenseMatrix & xg, unsigned int & gp){
    std::cout << "Not using this method in this application.\n";
    Epetra_SerialDenseVector f(3);
    return f;
}
Epetra_SerialDenseVector phaseFieldProblem::get_forcing(double & x1, double & x2, double & x3, unsigned int & e_lid, unsigned int & gp){
    std::cout << "Not using this method in this application.\n";
    Epetra_SerialDenseVector f(3);
    return f;
}

void phaseFieldProblem::setup_dirichlet_conditions(){

}

void phaseFieldProblem::apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){

}
