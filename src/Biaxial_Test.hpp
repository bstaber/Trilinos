#ifndef COMPRESSIBLE_MOONEY_TRANSVERSE_ISOTROPIC_HPP
#define COMPRESSIBLE_MOONEY_TRANSVERSE_ISOTROPIC_HPP

#include "tensor_calculus.hpp"
#include "hyperelasticity_setup_pp.hpp"

class Biaxial_Test : public hyperelasticity_setup
{
public:
    
    Biaxial_Test(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
        
        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh = new mesh(comm, mesh_file);
        Comm = Mesh->Comm;
        
        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();
        
        setup_dirichlet_conditions();
    }
    
    ~Biaxial_Test(){
    }
    
    void set_parameters(Epetra_SerialDenseVector & x){
        alpha1 = x(0);
        alpha2 = x(1);
        alpha3 = x(2);
        alpha4 = x(3);
        alpha5 = x(4);
        delta1 = x(5);
        eta1   = x(6);
        
        trm = mu4 + 2.0*mu5;
        ptrmbeta4 = std::pow(trm,alpha4);
        ptrmbeta5 = std::pow(trm,alpha5);
        
        beta1  = alpha4*ptrmbeta4/(alpha5*ptrmbeta5);
        delta2 = 2.0*alpha1 + 4.0*alpha2 + 2.0*delta1 + 2.0*eta1*alpha4*ptrmbeta4;
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
    
    void get_material_parameters(unsigned int & e_lid, unsigned int & gp){
    }
    
    void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola){
        
        double det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);
        
        double scalarAB = 4.0*alpha2;
        tensor_product(scalarAB,eye,eye,tangent_piola,0.0);
        scalarAB = -4.0*alpha2;
        sym_tensor_product(scalarAB,eye,eye,tangent_piola,1.0);
        scalarAB = 4.0*delta1*det*det;
        tensor_product(scalarAB,L,L,tangent_piola,1.0);
        scalarAB = 2.0*delta2 - 4.0*delta1*det*det;
        sym_tensor_product(scalarAB,L,L,tangent_piola,1.0);
        
        scalarAB = 4.0*eta1*alpha4*(alpha4-1.0)*pJ4/(J4*J4);
        tensor_product(scalarAB,M,M,tangent_piola,1.0);
        scalarAB = 4.0*eta1*beta1*alpha5*(alpha5-1.0)*pJ5/J5;
        tensor_product(scalarAB,dJ5,dJ5,tangent_piola,1.0);
        
        scalarAB = J5;
        tensor_product(scalarAB,L,L,ddJ5,1.0);
        scalarAB = -J5;
        sym_tensor_product(scalarAB,L,L,ddJ5,1.0);
        scalarAB = -I3;
        tensor_product(scalarAB,L,LML,ddJ5,1.0);
        tensor_product(scalarAB,LML,L,ddJ5,1.0);
        scalarAB = I3;
        sym_tensor_product(scalarAB,L,LML,ddJ5,1.0);
        sym_tensor_product(scalarAB,LML,L,ddJ5,1.0);
        
        
        
    }
    
    void get_constitutive_tensors_static_condensation(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseVector & inverse_cauchy, Epetra_SerialDenseVector & piola_isc, Epetra_SerialDenseVector & piola_vol, Epetra_SerialDenseMatrix & tangent_piola_isc, Epetra_SerialDenseMatrix & tangent_piola_vol){
        std::cerr << "**Err: Not using static condensation method!\n";
    }
    
    void get_internal_pressure(double & theta, double & pressure, double & dpressure){
        std::cerr << "**Err: Not using static condensation method!\n";
    }
    
    void get_material_parameters_for_recover(unsigned int & e_lid, double & xi, double & eta, double & zeta){
    }
                                                      
    void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress){
    }
    
    double alpha1;
    double alpha2;
    double alpha3;
    double alpha4;
    double alpha5;
    double delta1;
    double delta2;
    double eta1;
    
    double trm;
    double ptrmbeta4;
    double ptrmbeta5;
    
};

#endif
