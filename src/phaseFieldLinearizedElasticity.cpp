/*
Brian Staber (brian.staber@gmail.com)
*/

#include "phaseFieldLinearizedElasticity.hpp"

phaseFieldLinearizedElasticity::phaseFieldLinearizedElasticity(Epetra_Comm & comm, Teuchos::ParameterList & Parameters,
                                                               double & gc_, double & lc_): gc(gc_), lc(lc_){

  std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
  Mesh = new mesh(comm, mesh_file, 1.0);
  Comm = Mesh->Comm;

  damageInterface = Teuchos::rcp(new damageField(comm, *Mesh, gc, lc));

  StandardMap = new Epetra_Map(-1, 3*Mesh->n_local_nodes_without_ghosts, &Mesh->local_dof_without_ghosts[0], 0, *Comm);
  OverlapMap  = new Epetra_Map(-1, 3*Mesh->n_local_nodes,&Mesh->local_dof[0], 0, *Comm);
  ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);
  create_FECrsGraph();
  constructGaussMap();

  damageHistory = new Epetra_Vector(*GaussMap);
  displacement  = new Epetra_Vector(*StandardMap);

  matrix        = new Epetra_FECrsMatrix(Copy,*FEGraph);
  rhs           = new Epetra_FEVector(*StandardMap);

  damageHistory->PutScalar(0.0);
  displacement->PutScalar(0.0);

  elasticity.Reshape(6,6);
  double c11 = E*(1.0-nu)/((1.0+nu)*(1.0-2.0*nu));
  double c12 = E*nu/((1.0+nu)*(1.0-2.0*nu));
  double c44 = E/(2.0*(1.0+nu));

  elasticity(0,0) = c11; elasticity(0,1) = c12; elasticity(0,2) = c12; elasticity(0,3) = 0.0; elasticity(0,4) = 0.0; elasticity(0,5) = 0.0;
  elasticity(1,0) = c12; elasticity(1,1) = c11; elasticity(1,2) = c12; elasticity(1,3) = 0.0; elasticity(1,4) = 0.0; elasticity(1,5) = 0.0;
  elasticity(2,0) = c12; elasticity(2,1) = c12; elasticity(2,2) = c11; elasticity(2,3) = 0.0; elasticity(2,4) = 0.0; elasticity(2,5) = 0.0;
  elasticity(3,0) = 0.0; elasticity(3,1) = 0.0; elasticity(3,2) = 0.0; elasticity(3,3) = c44; elasticity(3,4) = 0.0; elasticity(3,5) = 0.0;
  elasticity(4,0) = 0.0; elasticity(4,1) = 0.0; elasticity(4,2) = 0.0; elasticity(4,3) = 0.0; elasticity(4,4) = c44; elasticity(4,5) = 0.0;
  elasticity(5,0) = 0.0; elasticity(5,1) = 0.0; elasticity(5,2) = 0.0; elasticity(5,3) = 0.0; elasticity(5,4) = 0.0; elasticity(5,5) = c44;
}

phaseFieldLinearizedElasticity::~phaseFieldLinearizedElasticity(){

}

void phaseFieldLinearizedElasticity::computeDisplacement(){

  matrix->PutScalar(0.0);
  rhs->PutScalar(0.0);

  assemblePureDirichlet_homogeneousForcing(*matrix);
  apply_dirichlet_conditions(*matrix, *rhs, bc_disp);

  displacement->PutScalar(0.0);

  Epetra_LinearProblem problem;
  AztecOO solver;

  problem.SetOperator(matrix);
  problem.SetLHS(displacement);
  problem.SetRHS(rhs);

  solver.SetProblem(problem);
  //solver.SetParameters(*Krylov);

  solver.Iterate(2000, 1.0e-6);
}

void phaseFieldLinearizedElasticity::updateDamageHistory(){

  Epetra_Vector u(*OverlapMap);
  u.Import(*displacement, *ImportToOverlapMap, Insert);

  int n_gauss_points = Mesh->n_gauss_cells;

  double trepsilon, trepsilon2, potential;

  Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
  Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
  Epetra_SerialDenseVector cells_u(3*Mesh->el_type);
  Epetra_SerialDenseVector epsilon(6);

  int egid;
  for (unsigned int elid=0; elid<Mesh->n_local_cells; ++elid){
    egid = Mesh->local_cells[elid];
    node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
    for (unsigned int ddl=0; ddl<3; ++ddl){
      cells_u(3*inode+dll) = u[OverlapMap->LID(3*node+ddl)];
    }
    for (unsigned int gp=0; gp<n_gauss_points; ++gp){
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,elid);
            dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,elid);
            dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,elid);
        }
        compute_B_matrices(dx_shape_functions,matrix_B);
        epsilon.Multiply('N','N',1.0,matrix_B,cells_u,0.0);
        trepsilon  = epsilon(0) + epsilon(1) + epsilon(2);
        trepsilon2 = epsilon(0)*epsilon(0) + epsilon(1)*epsilon(1) + epsilon(2)*epsilon(2) +
                     0.5*epsilon(3)*epsilon(3) + 0.5*epsilon(4)*epsilon(4) + 0.5*epsilon(5)*epsilon(5);
        potential = (lambda/2.0)*trepsilon*trepsilon + mu*trepsilon2;
        if (potential>damageHistory[GaussMap->LID(n_gauss_points*egid+gp)]){
          damageHistory[GaussMap->LID(n_gauss_points*egid+gp)] = potential;
        }
    }
  }

}

void phaseFieldLinearizedElasticity::get_elasticity_tensor(unsigned int & e_lid,
                                                           unsigned int & gp,
                                                           Epetra_SerialDenseMatrix & tangent_matrix){
  int n_gauss_points = Mesh->n_gauss_cells;
  int e_gid = Mesh->local_cells[e_lid];
  int node;
  double xi = Mesh->xi_cells[gp]; double eta = Mesh->eta_cells[gp]; double zeta = Mesh->zeta_cells[gp];

  Epetra_SerialDenseVector shape_functions(Mesh->el_type);
  switch (Mesh->el_type){
      case 4:
          tetra4::shape_functions(shape_functions, xi, eta, zeta);
          break;
      case 8:
          hexa8::shape_functions(shape_functions, xi, eta, zeta);
          break;
      case 10:
          tetra10::shape_functions(shape_functions, xi, eta, zeta);
          break;
  }

  double d = 0.0;
  for (unsigned int j=0; j<Mesh->el_type; ++j){
      node = Mesh->cells_nodes[Mesh->el_type*e_gid+j];
      d += shape_functions(j)*damageInterface->damageSolution[0][damageInterface->OverlapMap->LID(node)];
  }

  double g = (1.0-d)*(1.0-d) + 1.0e-6;
  tangent_matrix = elasticity;
  tangent_matrix.Scale(g);
}

void phaseFieldLinearizedElasticity::get_elasticity_tensor_for_recovery(unsigned int & e_lid,
                                                                        Epetra_SerialDenseMatrix & tangent_matrix){
  std::cout << "Not using this method yet.\n";
}

void phaseFieldLinearizedElasticity::constructGaussMap(){
  int e_gid;
  int n_local_cells = interface->Mesh->n_local_cells;
  int n_gauss_cells = interface->Mesh->n_gauss_cells;
  std::vector<int> local_gauss_points(n_local_cells*n_gauss_cells);
  for (unsigned int e_lid=0; e_lid<n_local_cells; ++e_lid){
      e_gid = interface->Mesh->local_cells[e_lid];
      for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
          local_gauss_points[e_lid*n_gauss_cells+gp] = e_gid*n_gauss_cells+gp;
      }

  }
  GaussMap = new Epetra_Map(-1, n_local_cells*n_gauss_cells, &local_gauss_points[0], 0, Comm);
}
