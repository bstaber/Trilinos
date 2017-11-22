#ifndef MANUFACTUREDSOLUTION_HPP
#define MANUFACTUREDSOLUTION_HPP

#include "tensor_calculus.hpp"
#include "Newton_Raphsonpp.hpp"
#include "distributenrldata.hpp"
#include "compressibleHyperelasticity.hpp"

class manufacturedSolution : public compressibleHyperelasticity
{
public:
    
    double mu1, mu2, mu3, mu4, mu5, mu;
    double trm, beta3, beta4, beta5;
    double ptrmbeta4, ptrmbeta5;
    double plyagl, cos_plyagl, sin_plyagl, topcoord, stepInc;
    
    manufacturedSolution(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
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
    
    ~manufacturedSolution(){
    }
    
    void get_matrix_and_rhs(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
        assembleMixedDirichletNeumann_inhomogeneousForcing(x,K,F);
    }
    
    void setup_dirichlet_conditions(){
        n_bc_dof = 0;
        double y;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            y = Mesh->nodes_coord[3*node+1];
            if(y==0.0){
                n_bc_dof+=3;
            }
            if(y==topcoord){
                n_bc_dof+=3;
            }
        }
        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            y = Mesh->nodes_coord[3*node+1];
            if (y==0.0){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
            if (y==topcoord){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
        }
    }
    
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        
        Epetra_MultiVector v(*StandardMap,true);
        Epetra_SerialDenseVector u(3);
        v.PutScalar(0.0);
        double x,y,z;
        int node;
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            if (y==0.0){
                u = getManufacturedSolution(x,y,z);
                v[0][StandardMap->LID(3*node+0)] = u(0);
                v[0][StandardMap->LID(3*node+1)] = u(1);
                v[0][StandardMap->LID(3*node+2)] = u(2);
            }
            if (y==topcoord){
                u = getManufacturedSolution(x,y,z);
                v[0][StandardMap->LID(3*node+0)] = u(0);
                v[0][StandardMap->LID(3*node+1)] = u(1);
                v[0][StandardMap->LID(3*node+2)] = u(2);
            }
        }
        Epetra_MultiVector rhs_dir(*StandardMap,true);
        K.Apply(v,rhs_dir);
        F.Update(-1.0,rhs_dir,1.0);
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            y = Mesh->nodes_coord[3*node+1];
            if (y==0.0){
                F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
                F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
                F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
            }
            if (y==topcoord){
                F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
                F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
                F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }
    
    Epetra_SerialDenseVector get_neumannBc(Epetra_SerialDenseMatrix & matrix_X, Epetra_SerialDenseMatrix & xg, unsigned int & gp){
        Epetra_SerialDenseVector h(3), normal(3);
        Epetra_SerialDenseMatrix piola(3,3), d_shape_functions(Mesh->face_type,2), dxi_matrix_x(3,2);
        piola = manufacturedStress(xg(0,gp),xg(1,gp),xg(2,gp));
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            d_shape_functions(inode,0) = Mesh->D1_N_tri(gp,inode);
            d_shape_functions(inode,1) = Mesh->D2_N_tri(gp,inode);
        }
        dxi_matrix_x.Multiply('N','N',1.0,matrix_X,d_shape_functions,0.0);
        normal(0) = dxi_matrix_x(1,0)*dxi_matrix_x(2,1) - dxi_matrix_x(2,0)*dxi_matrix_x(1,1);
        normal(1) = dxi_matrix_x(2,0)*dxi_matrix_x(0,1) - dxi_matrix_x(0,0)*dxi_matrix_x(2,1);
        normal(2) = dxi_matrix_x(0,0)*dxi_matrix_x(1,1) - dxi_matrix_x(1,0)*dxi_matrix_x(0,1);
        if (xg(2,gp)==0.0){
            normal.Scale(-1.0);
        }
        h.Multiply('N','N',1.0,piola,normal,0.0);
        return h;
    }
    
    Epetra_SerialDenseVector get_forcing(double & x1, double & x2, double & x3, unsigned int & e_lid, unsigned int & gp){
        Epetra_SerialDenseVector f(3);
        f = manufacturedForcing(x1,x2,x3);
        return f;
    }
    
    void get_material_parameters(unsigned int & e_lid, unsigned int & gp){
        //std::cout << "**Err: Not using that method in this example!\n";
    }
    
    void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola){
        
        double det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);
        
        Epetra_SerialDenseMatrix C(3,3);
        Epetra_SerialDenseMatrix CC(3,3);
        Epetra_SerialDenseMatrix ddJ5(6,6);
        Epetra_SerialDenseVector LML(6);
        Epetra_SerialDenseVector eye(6);
        Epetra_SerialDenseVector dJ5(6);
        Epetra_SerialDenseVector M(6);
        Epetra_SerialDenseVector L(6);
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
        
        M(0) = mu4*sin_plyagl*sin_plyagl+mu5*cos_plyagl*cos_plyagl;
        M(1) = mu4*cos_plyagl*cos_plyagl+mu5*sin_plyagl*sin_plyagl;
        M(2) = mu5;
        M(3) = 0.0;
        M(4) = 0.0;
        M(5) = (mu4-mu5)*cos_plyagl*sin_plyagl;
        
        LML(0) = M(0)*L(0)*L(0)+2.0*M(4)*L(0)*L(4)+2.0*M(5)*L(0)*L(5)+M(2)*L(4)*L(4)+2.0*M(3)*L(4)*L(5)+M(1)*L(5)*L(5);
        LML(1) = M(1)*L(1)*L(1)+2.0*M(4)*L(1)*L(3)+2.0*M(5)*L(1)*L(5)+M(2)*L(3)*L(3)+2.0*M(4)*L(3)*L(5)+M(0)*L(5)*L(5);
        LML(2) = M(2)*L(2)*L(2)+2.0*M(3)*L(2)*L(3)+2.0*M(4)*L(2)*L(4)+M(1)*L(3)*L(3)+2.0*M(5)*L(3)*L(4)+M(0)*L(4)*L(4);
        LML(3) = L(2)*(L(1)*M(3)+L(3)*M(2)+L(5)*M(4))+L(3)*(L(1)*M(1)+L(3)*M(3)+L(5)*M(5))+L(4)*(L(5)*M(0)+L(1)*M(5)+L(3)*M(4));
        LML(4) = L(2)*(L(0)*M(4)+L(4)*M(2)+L(5)*M(3))+L(3)*(L(0)*M(5)+L(5)*M(1)+L(4)*M(3))+L(4)*(L(0)*M(0)+L(4)*M(4)+L(5)*M(5));
        LML(5) = L(1)*(L(0)*M(5)+L(5)*M(1)+L(4)*M(3))+L(3)*(L(0)*M(4)+L(4)*M(2)+L(5)*M(3))+L(5)*(L(0)*M(0)+L(4)*M(4)+L(5)*M(5));
        
        eye(0) = 1.0; eye(1) = 1.0; eye(2) = 1.0; eye(3) = 0.0; eye(4) = 0.0; eye(5) = 0.0;
        
        double I1 = C(0,0) + C(1,1) + C(2,2);
        double II1 = C(0,0)*C(0,0) + C(1,1)*C(1,1) + C(2,2)*C(2,2) + 2.0*C(1,2)*C(1,2) + 2.0*C(0,2)*C(0,2) + 2.0*C(0,1)*C(0,1);
        double I2 = (1.0/2.0)*(I1*I1-II1);
        double I3 = det*det;
        double I4 = C(0,0)*M(0) + C(1,1)*M(1) + C(2,2)*M(2) + 2.0*C(0,1)*M(5) + 2.0*C(0,2)*M(4) + 2.0*C(1,2)*M(3);
        double I5 = CC(0,0)*M(0) + CC(1,1)*M(1) + CC(2,2)*M(2) + 2.0*CC(0,1)*M(5) + 2.0*CC(0,2)*M(4) + 2.0*CC(1,2)*M(3);
        double J5 = I5 - I1*I4 + I2*trm;
        double pI3 = std::pow(I3,-beta3);
        double pI4 = std::pow(I4,beta4);
        double pJ5 = std::pow(J5,beta5);
        
        for (unsigned int i=0; i<6; ++i){
            dJ5(i) = J5*L(i) - I3*LML(i);
            piola_stress(i) = 2.0*mu1*eye(i) + 2.0*mu2*(I1*eye(i)-c(i)) + (2.0*mu3*det*det-mu)*L(i) + (2.0/ptrmbeta4)*pI4*M(i) + (2.0/ptrmbeta5)*pJ5*dJ5(i) - 2.0*trm*pI3*L(i);
        }
        
        double scalarAB = 4.0*mu2;
        tensor_product(scalarAB,eye,eye,tangent_piola,0.0);
        
        scalarAB = -4.0*mu2;
        sym_tensor_product(scalarAB,eye,eye,tangent_piola,1.0);
        
        scalarAB = 4.0*mu3*det*det;
        tensor_product(scalarAB,L,L,tangent_piola,1.0);
        
        scalarAB = -4.0*mu3*det*det+2.0*mu;
        sym_tensor_product(scalarAB,L,L,tangent_piola,1.0);
        
        scalarAB = J5;
        tensor_product(scalarAB,L,L,ddJ5,0.0);
        scalarAB = -J5;
        sym_tensor_product(scalarAB,L,L,ddJ5,1.0);
        scalarAB = -I3;
        tensor_product(scalarAB,L,LML,ddJ5,1.0);
        tensor_product(scalarAB,LML,L,ddJ5,1.0);
        scalarAB = I3;
        sym_tensor_product(scalarAB,L,LML,ddJ5,1.0);
        sym_tensor_product(scalarAB,LML,L,ddJ5,1.0);
        
        scalarAB = 4.0*beta4*pI4/(I4*ptrmbeta4);
        tensor_product(scalarAB,M,M,tangent_piola,1.0);
        
        scalarAB = 4.0*beta5*pJ5/(J5*ptrmbeta5);
        tensor_product(scalarAB,dJ5,dJ5,tangent_piola,1.0);
        
        scalarAB = 4.0*trm*beta3*pI3;
        tensor_product(scalarAB,L,L,tangent_piola,1.0);
        scalarAB = 4.0*trm*pI3;
        sym_tensor_product(scalarAB,L,L,tangent_piola,1.0);
        
        ddJ5.Scale(4.0*pJ5/ptrmbeta5);
        tangent_piola += ddJ5;
    }
    
    
    Epetra_SerialDenseVector getManufacturedSolution(double & x1, double & x2, double & x3){
        Epetra_SerialDenseVector u(3);
        Epetra_SerialDenseVector u(3);
        double c1 = 2.0e-4; double c2 = 1.0e-4; double c3 = 2.0e-4; double topcoord = 25.0;
        u(0) = -c1*x2*(topcoord - x1)*(topcoord - x2);
        u(1) = c2*x2*(topcoord/2.0 - x2);
        u(2) = std::sin(c1*x3);
        return u;
    }
    
    Epetra_SerialDenseMatrix getManufacturedPiola(double & x1, double & x2, double & x3){
        
        Epetra_SerialDenseVector x(3);
        Epetra_SerialDenseMatrix F(3,3), C(3,3), CC(3,3), ML(3,3), LML(3,3), L(3,3), M(3,3), eye(3,3), S(3,3), P(3,3);
        double c1 = 2.0e2;
        double c2 = 1.0e1;
        double c3 = 2.0e2;
        x(0) = x1; x(1) = x2; x(2) = x3;
        F(0,0) = c1*x2*(topcoord-x2) + 1.0; F(0,1) = c1*x2*(topcoord-x1)-c1*(topcoord-x1)*(topcoord-x2); F(0,2) = 0.0;
        F(1,0) = 0.0;                       F(1,1) = c2*(topcoord/2.0 - x2)-c2*x2+1.0;                   F(1,2) = 0.0;
        F(2,0) = 0.0;                       F(2,1) = 0.0;                                                F(2,2) = c1*cos(c1*x3) + 1.0;
        double det = F(0,0)*F(1,1)*F(2,2)-F(0,0)*F(1,2)*F(2,1)-F(0,1)*F(1,0)*F(2,2)+F(0,1)*F(1,2)*F(2,0)+F(0,2)*F(1,0)*F(2,1)-F(0,2)*F(1,1)*F(2,0);
        
        M(0,0) = mu4*sin_plyagl*sin_plyagl+mu5*cos_plyagl*cos_plyagl;
        M(1,1) = mu4*cos_plyagl*cos_plyagl+mu5*sin_plyagl*sin_plyagl;
        M(2,2) = mu5;
        M(1,2) = 0.0; M(2,1) = 0.0;
        M(0,2) = 0.0; M(2,0) = 0.0;
        M(0,1) = (mu4-mu5)*cos_plyagl*sin_plyagl; M(1,0) = M(0,1);
        
        C.Multiply('T','N',1.0,F,F,0.0);
        CC.Multiply('N','N',1.0,C,C,0.0);
        
        L(0,0) = (1.0/(det*det))*(C(1,1)*C(2,2)-C(1,2)*C(2,1));
        L(1,1) = (1.0/(det*det))*(C(0,0)*C(2,2)-C(0,2)*C(2,0));
        L(2,2) = (1.0/(det*det))*(C(0,0)*C(1,1)-C(0,1)*C(1,0));
        L(1,2) = (1.0/(det*det))*(C(0,2)*C(1,0)-C(0,0)*C(1,2)); L(2,1) = L(1,2);
        L(0,2) = (1.0/(det*det))*(C(0,1)*C(1,2)-C(0,2)*C(1,1)); L(2,0) = L(0,2);
        L(0,1) = (1.0/(det*det))*(C(0,2)*C(2,1)-C(0,1)*C(2,2)); L(1,0) = L(0,1);
        
        ML.Multiply('N','N',1.0,M,L,0.0);
        LML.Multiply('N','N',1.0,L,ML,0.0);
        
        double I1 = C(0,0) + C(1,1) + C(2,2);
        double II1 = C(0,0)*C(0,0) + C(1,1)*C(1,1) + C(2,2)*C(2,2) + 2.0*C(1,2)*C(1,2) + 2.0*C(0,2)*C(0,2) + 2.0*C(0,1)*C(0,1);
        double I2 = (1.0/2.0)*(I1*I1-II1);
        double I3 = det*det;
        double I4 = C(0,0)*M(0,0) + C(1,1)*M(1,1) + C(2,2)*M(2,2) + 2.0*C(0,1)*M(0,1) + 2.0*C(0,2)*M(0,2) + 2.0*C(1,2)*M(1,2);
        double I5 = CC(0,0)*M(0,0) + CC(1,1)*M(1,1) + CC(2,2)*M(2,2) + 2.0*CC(0,1)*M(0,1) + 2.0*CC(0,2)*M(0,2) + 2.0*CC(1,2)*M(1,2);
        double J5 = I5 - I1*I4 + I2*trm;
        double pI3 = std::pow(I3,-beta3);
        double pI4 = std::pow(I4,beta4);
        double pJ5 = std::pow(J5,beta5);
        
        eye(0,0) = 1.0; eye(0,1) = 0.0; eye(0,2) = 0.0;
        eye(1,0) = 0.0; eye(1,1) = 1.0; eye(1,2) = 0.0;
        eye(2,0) = 0.0; eye(2,1) = 0.0; eye(2,2) = 1.0;
        
        for (unsigned int i=0; i<3; ++i){
            for (unsigned int j=0; j<3; ++j){
                S(i,j) = 2.0*mu1*eye(i,j) + 2.0*mu2*(I1*eye(i,j)-C(i,j)) + (2.0*mu3*det*det-mu)*L(i,j) + (2.0/ptrmbeta4)*pI4*M(i,j) + (2.0/ptrmbeta5)*pJ5*(J5*L(i,j)-I3*LML(i,j)) - 2.0*trm*pI3*L(i,j);
            }
        }
        P.Multiply('N','N',1.0,F,S,0.0);
        return P;
    }
    
    Epetra_SerialDenseVector manufacturedForcing(double & x1, double & x2, double & x3){
        Epetra_SerialDenseVector f(3), x(3), xf(3), xb(3);
        Epetra_SerialDenseMatrix Pf(3,3), Pb(3,3);
        double h = 1.0e-6;
        x(0) = x1; x(1) = x2; x(2) = x3;
        for (unsigned int j=0; j<3; ++j){
            xf = x;     xb = x;
            xf(j) += h; xb(j) -= h;
            Pf = getManufacturedPiola(xf(0),xf(1),xf(2));
            Pb = getManufacturedPiola(xb(0),xb(1),xb(2));
            for (unsigned int i=0; i<3; ++i){
                f(i) -= (Pf(i,j)-Pb(i,j))/(2.0*h);
            }
        }
        return f;
    }
    
    void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress){
        std::cout << "**Err: Not using that method in this example!\n";
    }
    
};

#endif

