#include "damageField.hpp"
#include "fepp.hpp"

damageField::damageField(Epetra_Comm & comm, mesh & mesh, double & gc_, double & lc_):
gc(gc_), lc(lc_){

  Mesh               = &mesh;
  Comm               = Mesh->Comm;

  StandardMap        = new Epetra_Map(-1, Mesh->n_local_nodes_without_ghosts, &Mesh->local_nodes_without_ghosts[0], 0, *Comm);
  OverlapMap         = new Epetra_Map(-1, Mesh->n_local_nodes, &Mesh->local_nodes[0], 0, *Comm);
  ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);

  create_FECrsGraph();

  damageSolution     = new Epetra_Vector(*StandardMap);
  matrix             = new Epetra_FECrsMatrix(Copy, *FEGraph);
  rhs                = new Epetra_FEVector(*StandardMap);
}

damageField::~damageField(){
}

void damageField::assemble(Epetra_Vector & Hn){

  //Hn to do

  Epetra_SerialDenseVector shape_functions(Mesh->el_type);
  Epetra_SerialDenseVector fe(Mesh->el_type);
  Epetra_SerialDenseVector hn(Mesh->el_type);

  Epetra_SerialDenseMatrix dx_shape_functions(3,Mesh->el_type);
  Epetra_SerialDenseMatrix ke(Mesh->el_type, Mesh->el_type);
  Epetra_SerialDenseMatrix me(Mesh->el_type, Mesh->el_type);

  double gauss_weight;
  unsigned int eglob;
  int n_gauss_points = Mesh->n_gauss_cells;
  int * index = new int [Mesh->el_type];

  for (unsigned int eloc=0; eloc<Mesh->n_local_cells; ++eloc){
    eglob = Mesh->local_cells[eloc];
    for (int inode=0; inode<Mesh->el_type; ++inode){
      index[inode] = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
      hn(inode) = Hn[OverlapMap->LID(inode)];
      fe(inode) = 0.0;
      for (int jnode=0; jnode<Mesh->el_type; ++jnode){
        ke(inode,jnode) = 0.0;
        me(inode,jnode) = 0.0;
      }
    }

    for (unsigned int gp=0; gp<n_gauss_points; ++gp){
      gauss_weight = Mesh->gauss_weight_cells(gp);
      double hn = 0.0; //todo
      double an = 0.0; //todo
      double bn = 0.0; //todo
      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          dx_shape_functions(0,inode) = Mesh->DX_N_cells(gp+n_gauss_points*inode,eloc);
          dx_shape_functions(1,inode) = Mesh->DY_N_cells(gp+n_gauss_points*inode,eloc);
          dx_shape_functions(2,inode) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,eloc);
          shape_functions(inode) = Mesh->N_cells(inode,gp);
          fe(inode) += gauss_weight*2.0*hn*shape_functions(inode)*Mesh->detJac_cells(eloc,gp);
      }
      me.Multiply('N','T',an*gauss_weight*Mesh->detJac_cells(eloc,gp),shape_functions,shape_functions,1.0);
      ke.Multiply('T','N',bn*gauss_weight*Mesh->detJac_cells(eloc,gp),dx_shape_functions,dx_shape_functions,1.0);
    }
    ke += me;

    for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
      rhs->SumIntoGlobalValues(1, &index[inode], &fe(inode));
      for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
        matrix->SumIntoGlobalValues(1, &index[inode], 1, &index[jnode], &ke(inode,jnode));
      }
    }
  }

  delete [] index;
}

void damageField::solve(Teuchos::ParameterList & Parameters,
                        Epetra_Vector & Hn){

  damageSolution->PutScalar(0.0);
  rhs->PutScalar(0.0);
  matrix->PutScalar(0.0);

  assemble(Hn);

  Epetra_LinearProblem problem(matrix, damageSolution, rhs);
  AztecOO solver(problem);

  double tol   = Teuchos::getParameter<double>(Parameters.sublist("Aztec"), "AZ_tol");
  int max_iter = Teuchos::getParameter<int>(Parameters.sublist("Aztec"), "AZ_max_iter");

  solver.SetParameters(Parameters);
  solver.Iterate(max_iter, tol);
}

void damageField::create_FECrsGraph(){
  FEGraph = new Epetra_FECrsGraph(Copy,*StandardMap,100);
  int eglob, node;
  int *index;
  index = new int [Mesh->el_type];

  for (int eloc=0; eloc<Mesh->n_local_cells; ++eloc){
      eglob = Mesh->local_cells[eloc];
      for (int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
          index[inode] = node;
      }

      for (int i=0; i<Mesh->el_type; ++i){
          for (int j=0; j<Mesh->el_type; ++j){
              FEGraph->InsertGlobalIndices(1, &index[i], 1, &index[j]);
          }
      }
  }
  Comm->Barrier();
  FEGraph->GlobalAssemble();
  delete[] index;
}

void damageField::setup_dirichlet_conditions(){
  std::cout << "No essential boundary conditions.\n";
}

void damageField::apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
  std::cout << "No essential boundary conditions.\n";
}
