#ifndef NRL_Neo_Hookean_hpp
#define NRL_Neo_Hookean_hpp

#include "tensor_calculus.hpp"
#include "hyperelasticity_setup_pp.hpp"

class Interface_dirichlet : public hyperelasticity_setup
{
public:
    
    Interface_dirichlet(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
        
        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh = new mesh(comm, mesh_file);
        Comm = Mesh->Comm;
        
        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();
        
        params.Reshape(2,1);
        
        setup_dirichlet_conditions();
    }
    
    ~Interface_dirichlet(){
    }
    
    void get_matrix_and_rhs(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
        assemble_dirichlet(x,K,F);
    }
    
    void setup_dirichlet_conditions(){
        
        n_bc_dof = 0;
        int dof = 1;
        double coord;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            coord = Mesh->nodes_coord[3*node+dof];
            if(coord==0){
                n_bc_dof+=3;
            }
            if(coord==25){
                n_bc_dof+=1;
            }
        }
        
        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+dof];
            if (coord==0){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
            if (coord==25){
                dof_on_boundary[indbc] = 3*inode+dof;
                indbc+=1;
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
                if (coord==25){
                    v[0][StandardMap->LID(3*node+dof)] = displacement;
                }
            }
        
            Epetra_MultiVector rhs_dir(*StandardMap,true);
            K.Apply(v,rhs_dir);
            F.Update(-1.0,rhs_dir,1.0);
        
            for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
                node = Mesh->local_nodes[inode];
                coord = Mesh->nodes_coord[3*node+dof];
                if (coord==0){
                    F[0][StandardMap->LID(3*node+0)] = 0.0;
                    F[0][StandardMap->LID(3*node+1)] = 0.0;
                    F[0][StandardMap->LID(3*node+2)] = 0.0;
                }
                if (coord==25){
                    F[0][StandardMap->LID(3*node+dof)] = displacement;
                }
            }
        //}
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }
    
    void get_material_parameters(unsigned int & e_lid, unsigned int & gp){
        params(0,0) = 10.0;
        params(1,0) = 50.0;
    }
    
    void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola){
        
        double det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);
        
        Epetra_SerialDenseMatrix rightCauchy(3,3);
        Epetra_SerialDenseVector eye(6);
        Epetra_SerialDenseVector L(6);
        
        rightCauchy.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);
        
        L(0) = (1.0/(det*det))*(rightCauchy(1,1)*rightCauchy(2,2)-rightCauchy(1,2)*rightCauchy(2,1));
        L(1) = (1.0/(det*det))*(rightCauchy(0,0)*rightCauchy(2,2)-rightCauchy(0,2)*rightCauchy(2,0));
        L(2) = (1.0/(det*det))*(rightCauchy(0,0)*rightCauchy(1,1)-rightCauchy(0,1)*rightCauchy(1,0));
        L(3) = (1.0/(det*det))*(rightCauchy(0,2)*rightCauchy(1,0)-rightCauchy(0,0)*rightCauchy(1,2));
        L(4) = (1.0/(det*det))*(rightCauchy(0,1)*rightCauchy(1,2)-rightCauchy(0,2)*rightCauchy(1,1));
        L(5) = (1.0/(det*det))*(rightCauchy(0,2)*rightCauchy(2,1)-rightCauchy(0,1)*rightCauchy(2,2));
        
        eye(0) = 1.0; eye(1) = 1.0; eye(2) = 1.0; eye(3) = 0.0; eye(4) = 0.0; eye(5) = 0.0;
        
        for (unsigned int i=0; i<6; ++i){
            piola_stress(i) = params(1,0)*eye(i) + (params(0,0)*det*(det-1.0)-params(1,0))*L(i);
        }
        
        double scalarAB = params(0,0)*det*(2.0*det-1.0);
        tensor_product(scalarAB,L,L,tangent_piola,0.0);
        scalarAB = 2.0*params(1,0)-params(0,0)*det*(det-1.0);
        sym_tensor_product(scalarAB,L,L,tangent_piola,1.0);
    }
    
    void get_constitutive_tensors_static_condensation(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseVector & piola_isc, Epetra_SerialDenseVector & piola_vol, Epetra_SerialDenseMatrix & tangent_piola_isc, Epetra_SerialDenseMatrix & tangent_piola_vol){
        std::cerr << "**Err: Not using static condensation method!\n";
    }
    
    void get_internal_pressure(double & theta, double & pressure, double & dpressure){
        std::cerr << "**Err: Not using static condensation method!\n";
    }
    
    void get_material_parameters_for_recover(unsigned int & e_lid, double & xi, double & eta, double & zeta){
    }
    void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress){
    }
    
    Epetra_SerialDenseMatrix params;
    
};

#endif
