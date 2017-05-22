#ifndef NRL_ModelB_hpp
#define NRL_ModelB_hpp

#include "tensor_calculus.hpp"
#include "hyperelasticity_setup_pp.hpp"

class Interface_dirichlet : public hyperelasticity_setup
{
public:
    
    Interface_dirichlet(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
        
        mu1 = Teuchos::getParameter<double>(Parameters.sublist("Model"), "mu1");
        mu4 = Teuchos::getParameter<double>(Parameters.sublist("Model"), "mu4");
        beta4 = Teuchos::getParameter<double>(Parameters.sublist("Model"), "beta4");
        beta5 = Teuchos::getParameter<double>(Parameters.sublist("Model"), "beta5");
        epsilon1 = Teuchos::getParameter<double>(Parameters.sublist("Model"), "epsilon1");
        
        mu5 = mu4*beta4/beta5;
        epsilon2 = 2.0*mu1 + 2.0*epsilon1 + 2.0*mu5*beta5;
        
        plyagl = Teuchos::getParameter<double>(Parameters.sublist("Model"), "angle");
        cos_plyagl = cos(plyagl);
        sin_plyagl = sin(plyagl);
        
        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        n_ply = Teuchos::getParameter<int>(Parameters.sublist("Mesh"), "n_ply");
        Mesh = new mesh(comm, mesh_file);
        Comm = Mesh->Comm;
        
        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();
        
        setup_dirichlet_conditions();
        for (unsigned int e=0; e<Mesh->n_cells/n_ply; ++e){
            for (unsigned int j=0; j<n_ply/2; ++j){
                phase.push_back(0);
                phase.push_back(1);
            }
        }
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
    }
    
    void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola){
        
        double det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);
        
        double pressure = 2.0*epsilon1*det - epsilon2/det;
        double dpressure = 2.0*epsilon1 + epsilon2/det;
        
        Epetra_SerialDenseMatrix C(3,3);
        Epetra_SerialDenseMatrix CC(3,3);
        Epetra_SerialDenseMatrix ddJ5(6,6);
        Epetra_SerialDenseVector LML(6);
        Epetra_SerialDenseVector eye(6);
        Epetra_SerialDenseVector dJ5(6);
        Epetra_SerialDenseVector M(6);
        Epetra_SerialDenseVector L(6);
        Epetra_SerialDenseVector D(6);
        Epetra_SerialDenseVector c(6);
        
        C.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);
        CC.Multiply('N','N',1.0,C,C,0.0);
        
        c(0) = C(0,0); c(1) = C(1,1); c(2) = C(2,2); c(3) = C(1,2); c(4) = C(0,2); c(5) = C(0,1);
        
        L(0) = (1.0/(det*det))*(C(1,1)*C(2,2)-C(1,2)*C(2,1));
        L(1) = (1.0/(det*det))*(C(0,0)*C(2,2)-C(0,2)*C(2,0));
        L(2) = (1.0/(det*det))*(C(0,0)*C(1,1)-C(0,1)*C(1,0));
        L(3) = (1.0/(det*det))*(C(0,2)*C(1,0)-C(0,0)*C(1,2));
        L(4) = (1.0/(det*det))*(C(0,1)*C(1,2)-C(0,2)*C(1,1));
        L(5) = (1.0/(det*det))*(C(0,2)*C(2,1)-C(0,1)*C(2,2));
        
        M(0) = cos_plyagl*cos_plyagl; M(1) = sin_plyagl*sin_plyagl; M(2) = 0.0; M(3) = 0.0; M(4)= 0.0; M(5) = cos_plyagl*sin_plyagl;
        
        LML(0) = M(0)*L(0)*L(0) + 2.0*M(4)*L(0)*L(4) + 2.0*M(5)*L(0)*L(5) + M(2)*L(4)*L(4) + 2.0*M(3)*L(4)*L(5) + M(1)*L(5)*L(5);
        LML(1) = M(1)*L(1)*L(1) + 2.0*M(4)*L(1)*L(3) + 2.0*M(5)*L(1)*L(5) + M(2)*L(3)*L(3) + 2*M(4)*L(3)*L(5) + M(0)*L(5)*L(5);
        LML(2) = M(2)*L(2)*L(2) + 2.0*M(3)*L(2)*L(3) + 2.0*M(4)*L(2)*L(4) + M(1)*L(3)*L(3) + 2.0*M(5)*L(3)*L(4) + M(0)*L(4)*L(4);
        LML(3) = L(2)*(L(1)*M(3) + L(3)*M(2) + L(5)*M(4)) + L(3)*(L(1)*M(1) + L(3)*M(3) + L(5)*M(5)) + L(4)*(L(5)*M(0) + L(1)*M(5) + L(3)*M(4));
        LML(4) = L(2)*(L(0)*M(4) + L(4)*M(2) + L(5)*M(3)) + L(3)*(L(0)*M(5) + L(5)*M(1) + L(4)*M(3)) + L(4)*(L(0)*M(0) + L(4)*M(4) + L(5)*M(5));
        LML(5) = L(1)*(L(0)*M(5) + L(5)*M(1) + L(4)*M(3)) + L(3)*(L(0)*M(4) + L(4)*M(2) + L(5)*M(3)) + L(5)*(L(0)*M(0) + L(4)*M(4) + L(5)*M(5));
        
        eye(0) = 1.0; eye(1) = 1.0; eye(2) = 1.0; eye(3) = 0.0; eye(4) = 0.0; eye(5) = 0.0;
        
        double I1 = C(0,0) + C(1,1) + C(2,2);
        double II1 = C(0,0)*C(0,0) + C(1,1)*C(1,1) + C(2,2)*C(2,2) + 2.0*C(1,2)*C(1,2) + 2.0*C(0,2)*C(0,2) + 2.0*C(0,1)*C(0,1);
        double I2 = (1.0/2.0)*(I1*I1-II1);
        double I3 = det*det;
        double I4 = C(0,0)*M(0) + C(1,1)*M(1) + C(2,2)*M(2) + 2.0*C(0,1)*M(5) + 2.0*C(0,2)*M(4) + 2.0*C(1,2)*M(3);
        double I5 = CC(0,0)*M(0) + CC(1,1)*M(1) + CC(2,2)*M(2) + 2.0*CC(0,1)*M(5) + 2.0*CC(0,2)*M(4) + 2.0*CC(1,2)*M(3);
        double J5 = I5 - I1*I4 + I2;
        double pI4 = std::pow(I4,beta4-1.0);
        double pJ5 = std::pow(J5,beta5-1.0);
        
        for (unsigned int i=0; i<6; ++i){
            dJ5(i) = J5*L(i) - I3*LML(i);
            piola_stress(i) = 2.0*mu1*eye(i) + det*pressure*L(i) + 2.0*mu4*beta4*pI4*M(i) + 2.0*mu5*beta5*pJ5*dJ5(i);
        }
        
        double scalarAB = J5;
        tensor_product(scalarAB,L,L,ddJ5,0.0);
        scalarAB = -J5;
        sym_tensor_product(scalarAB,L,L,ddJ5,1.0);
        scalarAB = -I3;
        tensor_product(scalarAB,L,LML,ddJ5,1.0);
        tensor_product(scalarAB,LML,L,ddJ5,1.0);
        scalarAB = I3;
        sym_tensor_product(scalarAB,L,LML,ddJ5,1.0);
        sym_tensor_product(scalarAB,LML,L,ddJ5,1.0);
        
        scalarAB = det*(pressure + det*dpressure);
        tensor_product(scalarAB,L,L,tangent_piola,0.0);
        scalarAB = -2.0*det*pressure;
        sym_tensor_product(scalarAB,L,L,tangent_piola,1.0);
        
        scalarAB = 4.0*mu4*beta4*(beta4-1.0)*(pI4/I4);
        tensor_product(scalarAB,M,M,tangent_piola,1.0);
        scalarAB = 4.0*mu5*beta5*(beta5-1.0)*(pJ5/J5);
        tensor_product(scalarAB,dJ5,dJ5,tangent_piola,1.0);
        ddJ5.Scale(4.0*mu5*beta5*pJ5);
        tangent_piola += ddJ5;
        
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
    
    double mu1;
    double epsilon1;
    double epsilon2;
    double mu4;
    double mu5;
    
    double beta4;
    double beta5;
    
    double plyagl;
    double cos_plyagl;
    double sin_plyagl;
    
    int n_ply;
    std::vector<int> phase;
};

#endif
