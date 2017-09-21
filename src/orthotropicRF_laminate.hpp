#ifndef ORTHOTROPICRF_LAMINATE_HPP
#define ORTHOTROPICRF_LAMINATE_HPP

#include "shinozukapp.hpp"
#include "tensor_calculus.hpp"
#include "linearelasticity_setup_pp.hpp"

class OrthotropicRF_Laminate : public LinearizedElasticity
{
public:
    
    OrthotropicRF_Laminate(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
        
        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh = new mesh(comm, mesh_file);
        Comm = Mesh->Comm;
        
        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();
        
        setup_dirichlet_conditions();
    }
    
    ~OrthotropicRF_Laminate(){
    }
        
    void setup_dirichlet_conditions(){
        n_bc_dof = 0;
        int dof = 1;
        double coord;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            coord = Mesh->nodes_coord[3*node+dof];
            if(coord==0.0){
                n_bc_dof+=3;
            }
            if(coord==25.0/1000.0){
                n_bc_dof+=3;
            }
        }
        
        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+dof];
            if (coord==0.0){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
            if (coord==25.0/1000.0){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+dof;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
        }
    }
    
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        //if (n_bc_dof>0){
        Epetra_MultiVector v(*StandardMap,true);
        v.PutScalar(0.0);
        
        int node;
        int dof = 1;
        double coord;
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+dof];
            if (coord==25.0/1000.0){
                v[0][StandardMap->LID(3*node+dof)] = displacement;
            }
        }
        
        Epetra_MultiVector rhs_dir(*StandardMap,true);
        K.Apply(v,rhs_dir);
        F.Update(-1.0,rhs_dir,1.0);
        
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+dof];
            if (coord==0.0){
                F[0][StandardMap->LID(3*node+0)] = 0.0;
                F[0][StandardMap->LID(3*node+1)] = 0.0;
                F[0][StandardMap->LID(3*node+2)] = 0.0;
            }
            if (coord==25.0/1000.0){
                F[0][StandardMap->LID(3*node+0)]   = 0.0;
                F[0][StandardMap->LID(3*node+dof)] = displacement;
                F[0][StandardMap->LID(3*node+2)]   = 0.0;
            }
        }
        //}
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }
    
    void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix){
        
        for (unsigned int i=0; i<tangent_matrix.M(); ++i){
            for (unsigned int j=0; j<tangent_matrix.N(); ++j){
                if (i==j){
                    tangent_matrix(i,j) = 1.0;
                }
                else{
                    tangent_matrix(i,j) = 0.0;
                }
            }
        }
        
    }
    
};

#endif
