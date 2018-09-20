#include "damageField.hpp"
#include "fepp.hpp"

damageField::damageField(mesh & mesh){
  Mesh = &mesh;
  Comm = Mesh->Comm;
  StandardMap        = new Epetra_Map(-1,Mesh->n_local_nodes_without_ghosts,&Mesh->local_nodes_without_ghosts[0],0,*Comm);
  OverlapMap         = new Epetra_Map(-1,Mesh->n_local_nodes,&Mesh->local_nodes[0],0,*Comm);
  ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);
  create_FECrsGraph();
}

damageField::~damageField(){
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
