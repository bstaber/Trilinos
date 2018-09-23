/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef PHASEFIELDPROBLEM_HPP
#define PHASEFIELDPROBLEM_HPP

#include "phaseFieldLinearizedElasticity.hpp"

class phaseFieldProblem : public phaseFieldLinearizedElasticity{

public:
  phaseFieldProblem(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
    initialize(comm, Parameters);
    setup_dirichlet_conditions();
  }
  ~phaseFieldProblem(){
  }

  Epetra_SerialDenseVector get_neumannBc(Epetra_SerialDenseMatrix & matrix_X, Epetra_SerialDenseMatrix & xg, unsigned int & gp){
      std::cout << "Not using this method in this application.\n";
      Epetra_SerialDenseVector f(3);
      return f;
  }
  Epetra_SerialDenseVector get_forcing(double & x1, double & x2, double & x3, unsigned int & e_lid, unsigned int & gp){
      std::cout << "Not using this method in this application.\n";
      Epetra_SerialDenseVector f(3);
      return f;
  }

  void setup_dirichlet_conditions(){

    n_bc_dof = 0;
    double y;
    int node;
    for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
        node = Mesh->local_nodes[i];
        y    = Mesh->nodes_coord[3*node+1];
        if(y==0.0){
            n_bc_dof+=3;
        }
        if(y==25.0){
            n_bc_dof+=1;
        }
    }

    int indbc = 0;
    dof_on_boundary = new int [n_bc_dof];
    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        y    = Mesh->nodes_coord[3*node+1];
        if (y==0.0){
            dof_on_boundary[indbc+0] = 3*inode+0;
            dof_on_boundary[indbc+1] = 3*inode+1;
            dof_on_boundary[indbc+2] = 3*inode+2;
            indbc+=3;
        }
        if (y==25.0){
            dof_on_boundary[indbc+0] = 3*inode+1;
            indbc+=1;
        }
    }

  }

  void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){

    Epetra_MultiVector v(*StandardMap,true);
    v.PutScalar(0.0);

    int node;
    double y;
    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        y    = Mesh->nodes_coord[3*node+1];
        if (y==25.0){
            v[0][StandardMap->LID(3*node+1)] = displacement;
        }
    }

    Epetra_MultiVector rhs_dir(*StandardMap,true);
    K.Apply(v,rhs_dir);
    F.Update(-1.0,rhs_dir,1.0);

    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        y    = Mesh->nodes_coord[3*node+1];
        if (y==0.0){
            F[0][StandardMap->LID(3*node+0)] = 0.0;
            F[0][StandardMap->LID(3*node+1)] = 0.0;
            F[0][StandardMap->LID(3*node+2)] = 0.0;
        }
        if (y==25.0){
            F[0][StandardMap->LID(3*node+1)] = displacement;
        }
    }
    //}
    ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
  }

};
#endif
