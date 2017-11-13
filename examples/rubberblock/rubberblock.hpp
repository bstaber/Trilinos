#ifndef COMPRESSIBLE_MOONEY_TRANSVERSE_ISOTROPIC_HPP
#define COMPRESSIBLE_MOONEY_TRANSVERSE_ISOTROPIC_HPP

#include "tensor_calculus.hpp"
#include "hyperelasticity_setup_pp.hpp"

class rubberblock : public hyperelasticity_setup
{
public:
    
    double topcoord;
    double lambda, mu;
    
    rubberblock(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
        
        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh = new mesh(comm, mesh_file);
        Comm = Mesh->Comm;
        findtop();
        
        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();
        
        setup_dirichlet_conditions();
    }
    
    ~rubberblock(){
    }
    
    void set_parameters(Epetra_SerialDenseVector & x){
        lambda = x(0);
        mu     = x(1);
    }
    
    void findtop(){
        topcoord = 0.0;
        if (Comm->MyPID()==0){
            for (unsigned int n=0; n<Mesh->n_nodes; ++n){
                if (Mesh->nodes_coord[3*n+1]>topcoord){
                    topcoord = Mesh->nodes_coord[3*n+1];
                }
            }
        }
        Comm->Broadcast(&topcoord,1,0);
        Comm->Barrier();
    }
        
    void get_matrix_and_rhs(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
        assemble_dirichlet(x,K,F);
    }
    
    void setup_dirichlet_conditions(){
        n_bc_dof = 0;
        double y;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            y    = Mesh->nodes_coord[3*node+1];
            if(y==0.0){
                n_bc_dof+=3;
            }
            if(y==topcoord){
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
            if (y==topcoord){
                dof_on_boundary[indbc] = 3*inode+1;
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
            y  = Mesh->nodes_coord[3*node+1];
            if (y==topcoord){
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
            if (y==topcoord){
                F[0][StandardMap->LID(3*node+1)] = displacement;
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }
    
    void get_material_parameters(unsigned int & e_lid, unsigned int & gp){
        //std::cerr << "**Err: Not using that method in this example!\n";
    }
    
    void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola){
        
        double det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);
        
        Epetra_SerialDenseMatrix C(3,3);
        Epetra_SerialDenseVector eye(6), L(6), c(6);
        
        C.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);
        c(0) = C(0,0); c(1) = C(1,1); c(2) = C(2,2); c(3) = C(1,2); c(4) = C(0,2); c(5) = C(0,1);
        
        L(0) = (1.0/(det*det))*(C(1,1)*C(2,2)-C(1,2)*C(2,1));
        L(1) = (1.0/(det*det))*(C(0,0)*C(2,2)-C(0,2)*C(2,0));
        L(2) = (1.0/(det*det))*(C(0,0)*C(1,1)-C(0,1)*C(1,0));
        L(3) = (1.0/(det*det))*(C(0,2)*C(1,0)-C(0,0)*C(1,2));
        L(4) = (1.0/(det*det))*(C(0,1)*C(1,2)-C(0,2)*C(1,1));
        L(5) = (1.0/(det*det))*(C(0,2)*C(2,1)-C(0,1)*C(2,2));
        
        eye(0) = 1.0; eye(1) = 1.0; eye(2) = 1.0; eye(3) = 0.0; eye(4) = 0.0; eye(5) = 0.0;
        
        double I1         = C(0,0) + C(1,1) + C(2,2);
        double dpressure  = (lambda/2.0)*det - ((lambda/2.0) + mu)/det;
        double ddpressure = (lambda/2.0)     + ((lambda/2.0) + mu)/(det*det);
        
        for (unsigned int i=0; i<6; ++i){
            piola_stress(i) = mu*eye(i) + det*dpressure*L(i);
        }
        
        double scalarAB = det*(dpressure+det*ddpressure);
        tensor_product(scalarAB,L,L,tangent_piola,0.0);
        scalarAB = -2.0*det*dpressure;
        sym_tensor_product(scalarAB,L,L,tangent_piola,1.0);
    }
    
    void compute_cauchy(Epetra_Vector & solution, double & xi, double & eta, double & zeta, std::string & filename){
        
        Epetra_Vector u(*OverlapMap);
        u.Import(solution, *ImportToOverlapMap, Insert);
        
        Epetra_Map CellsMap(-1,Mesh->n_local_cells,&Mesh->local_cells[0],0,*Comm);
        Epetra_Vector sig22(CellsMap);
        
        int node, e_gid;
        double det_jac_tetra;
        double I1, det, dpressure;
        
        Epetra_SerialDenseMatrix deformation_gradient(3,3), right_cauchy(3,3), inv_right_cauchy(3,3);
        Epetra_SerialDenseMatrix s(3,3), sig(3,3), sf(3,3), eye(3,3);
        Epetra_SerialDenseMatrix matrix_X(3,Mesh->el_type);
        Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
        Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
        Epetra_SerialDenseMatrix D(Mesh->el_type,3), DX(Mesh->el_type,3);
        Epetra_SerialDenseMatrix JacobianMatrix(3,3);
        
        eye(0,0) = 1.0; eye(0,1) = 0.0; eye(0,2) = 0.0;
        eye(1,0) = 0.0; eye(1,1) = 1.0; eye(1,2) = 0.0;
        eye(2,0) = 0.0; eye(2,1) = 0.0; eye(2,2) = 1.0;
        
        for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
            e_gid = Mesh->local_cells[e_lid];
            
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
                matrix_X(0,inode) = Mesh->nodes_coord[3*node+0];
                matrix_X(1,inode) = Mesh->nodes_coord[3*node+1];
                matrix_X(2,inode) = Mesh->nodes_coord[3*node+2];
                matrix_x(0,inode) = u[OverlapMap->LID(3*node+0)] + Mesh->nodes_coord[3*node+0];
                matrix_x(1,inode) = u[OverlapMap->LID(3*node+1)] + Mesh->nodes_coord[3*node+1];
                matrix_x(2,inode) = u[OverlapMap->LID(3*node+2)] + Mesh->nodes_coord[3*node+2];
            }
            
            switch (Mesh->el_type){
                case 4:
                    tetra4::d_shape_functions(D, xi, eta, zeta);
                    break;
                case 8:
                    hexa8::d_shape_functions(D, xi, eta, zeta);
                    break;
                case 10:
                    tetra10::d_shape_functions(D, xi, eta, zeta);
                    break;
            }
            jacobian_matrix(matrix_X,D,JacobianMatrix);
            jacobian_det(JacobianMatrix,det_jac_tetra);
            dX_shape_functions(D,JacobianMatrix,det_jac_tetra,dx_shape_functions);
            
            deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
            det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);
            right_cauchy.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);
            inv_right_cauchy(0,0) = (1.0/(det*det))*(right_cauchy(1,1)*right_cauchy(2,2)-right_cauchy(1,2)*right_cauchy(2,1));
            inv_right_cauchy(1,1) = (1.0/(det*det))*(right_cauchy(0,0)*right_cauchy(2,2)-right_cauchy(0,2)*right_cauchy(2,0));
            inv_right_cauchy(2,2) = (1.0/(det*det))*(right_cauchy(0,0)*right_cauchy(1,1)-right_cauchy(0,1)*right_cauchy(1,0));
            inv_right_cauchy(1,2) = (1.0/(det*det))*(right_cauchy(0,2)*right_cauchy(1,0)-right_cauchy(0,0)*right_cauchy(1,2));
            inv_right_cauchy(2,1) = inv_right_cauchy(1,2);
            inv_right_cauchy(0,2) = (1.0/(det*det))*(right_cauchy(0,1)*right_cauchy(1,2)-right_cauchy(0,2)*right_cauchy(1,1));
            inv_right_cauchy(2,0) = inv_right_cauchy(0,2);
            inv_right_cauchy(0,1) = (1.0/(det*det))*(right_cauchy(0,2)*right_cauchy(2,1)-right_cauchy(0,1)*right_cauchy(2,2));
            inv_right_cauchy(1,0) = inv_right_cauchy(0,1);
            I1         = right_cauchy(0,0) + right_cauchy(1,1) + right_cauchy(2,2);
            dpressure  = (lambda/2.0)*det - ((lambda/2.0) + mu)/det;
            for (unsigned int i=0; i<6; ++i){
                for (unsigned int j=0; j<6; ++j){
                    s(i,j) = mu*eye(i,j) + det*dpressure*inv_right_cauchy(i,j);
                }
            }
            sf.Multiply('N','T',1.0,s,deformation_gradient,0.0);
            sig.Multiply('N','N',1.0,deformation_gradient,sf,0.0);
            sig.Scale(1.0/det);
            //sig22[e_lid] = sig(1,1);
        }
        
        int NumTargetElements = 0;
        if (Comm->MyPID()==0){
            NumTargetElements = Mesh->n_cells;
        }
        Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
        Epetra_Export ExportOnRoot(CellsMap,MapOnRoot);
        Epetra_MultiVector lhs_root(MapOnRoot,true);
        lhs_root.Export(sig22,ExportOnRoot,Insert);
        
        int error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    }
    
    void get_constitutive_tensors_static_condensation(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseVector & inverse_cauchy, Epetra_SerialDenseVector & piola_isc, Epetra_SerialDenseVector & piola_vol, Epetra_SerialDenseMatrix & tangent_piola_isc, Epetra_SerialDenseMatrix & tangent_piola_vol){
        std::cerr << "**Err: Not using static condensation method!\n";
    }
    
    void get_internal_pressure(double & theta, double & pressure, double & dpressure){
        std::cerr << "**Err: Not using static condensation method!\n";
    }
    
    void get_material_parameters_for_recover(unsigned int & e_lid, double & xi, double & eta, double & zeta){
        std::cerr << "**Err: Not using that method in this example!\n";
    }
                                                      
    void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress){
        std::cerr << "**Err: Not using that method in this example!\n";
    }
    
};

#endif
