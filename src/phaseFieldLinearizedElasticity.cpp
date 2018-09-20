#include "phaseFieldLinearizedElasticity.hpp"

phaseFieldLinearizedElasticity::phaseFieldLinearizedElasticity(mesh & mesh, double & gc_, double & lc_): gc(gc_), lc(lc_){

  Mesh = &mesh;
  Comm = Mesh->Comm;

  damage = Teuchos::rcp(new damageField(mesh, gc, lc));

  StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
  OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
  ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
  create_FECrsGraph();

}

~phaseFieldLinearizedElasticity::phaseFieldLinearizedElasticity(){

}

void phaseFieldLinearizedElasticity::sovle(){

}

Epetra_SerialDenseVector phaseFieldLinearizedElasticity::get_neumannBc(Epetra_SerialDenseMatrix & matrix_X,
                                                                       Epetra_SerialDenseMatrix & xg,
                                                                       unsigned int & gp){

}

Epetra_SerialDenseVector phaseFieldLinearizedElasticity::get_forcing(double & x1, double & x2, double & x3,
                                                                     unsigned int & e_lid,
                                                                     unsigned int & gp){

}

void phaseFieldLinearizedElasticity::get_elasticity_tensor(unsigned int & e_lid,
                                                           unsigned int & gp,
                                                           Epetra_SerialDenseMatrix & tangent_matrix){

}

void phaseFieldLinearizedElasticity::get_elasticity_tensor_for_recovery(unsigned int & e_lid,
                                                                        Epetra_SerialDenseMatrix & tangent_matrix){

}

void phaseFieldLinearizedElasticity::setup_dirichlet_conditions(){

}

void phaseFieldLinearizedElasticity::apply_dirichlet_conditions(Epetra_FECrsMatrix & K,
                                                                Epetra_FEVector & F,
                                                                double & displacement){

}
