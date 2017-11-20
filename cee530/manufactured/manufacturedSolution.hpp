#ifndef manufactured_HPP
#define manufactured_HPP

#include "tensor_calculus.hpp"
#include "linearelasticity_setup_pp.hpp"

class manufactured : public LinearizedElasticity
{
public:
    
    Teuchos::ParameterList * Krylov;
    Epetra_Vector * uh;
    
    manufactured(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
        
        Krylov = &Parameters.sublist("Krylov");
        
        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh = new mesh(comm, mesh_file);
        Comm = Mesh->Comm;
        
        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();
        
        setup_dirichlet_conditions();
    }
    
    ~manufactured(){
    }
    
    Epetra_SerialDenseVector manufacturedSolution(double & x1, double & x2, double & x3){
        Epetra_SerialDenseVector u(3);
        u(0) = x1*x1*x1*x1 + (2.0*x2*x3)/5.0;
        u(1) = x2*x2*x2*x3 + (2.0*x1*x3)/5.0;
        u(2) = x1*x1;
        return u;
    }
    
    Epetra_SerialDenseMatrix manufacturedDeformation(double & x1, double & x2, double & x3){
        Epetra_SerialDenseMatrix epsilon(3,3);
        epsilon(0,0) = 4.0*x1*x1*x1;
        epsilon(1,1) = 4.0*x2*x2*x2;
        epsilon(2,2) = 0.0;
        epsilon(1,2) = x1/5.0; epsilon(2,1) = epsilon(1,2);
        epsilon(0,2) = x1 + x2/5.0; epsilon(2,0) = epsilon(0,2);
        epsilon(0,1) = (2.0*x3)/5.0; epsilon(1,0) = epsilon(0,1);
        return epsilon;
    }
    
    Epetra_SerialDenseMatrix manufacturedStress(double & x1, double & x2, double & x3){
        Epetra_SerialDenseMatrix sigma(3,3), epsilon(3,3), elasticity(6,6);
        Epetra_SerialDenseVector sigma_voigt(6), epsilon_voigt(6);
        epsilon = manufacturedDeformation(x1,x2,x3);
        epsilon_voigt(0) = epsilon(0,0); epsilon_voigt(1) = epsilon(1,1); epsilon_voigt(2) = epsilon(2,2);
        epsilon_voigt(3) = std::sqrt(2.0)*epsilon(1,2); epsilon_voigt(4) = std::sqrt(2.0)*epsilon(0,2); epsilon_voigt(5) = std::sqrt(2.0)*epsilon(0,1);
        unsigned int e_lid = 0;
        unsigned int gp    = 0;
        get_elasticity_tensor(e_lid,gp,elasticity);
        sigma_voigt.Multiply('N','N',1.0,elasticity,epsilon_voigt,0.0);
        sigma(0,0) = sigma_voigt(0);                sigma(0,1) = sigma_voigt(5)/std::sqrt(2.0); sigma(0,2) = sigma_voigt(4)/std::sqrt(2.0);
        sigma(1,0) = sigma_voigt(5)/std::sqrt(2.0); sigma(1,1) = sigma_voigt(1);                sigma(1,2) = sigma_voigt(3)/std::sqrt(2.0);
        sigma(2,0) = sigma_voigt(4)/std::sqrt(2.0); sigma(2,1) = sigma_voigt(3)/std::sqrt(2.0); sigma(2,2) = sigma_voigt(2);
        return sigma;
    }
    
    Epetra_SerialDenseVector get_neumannBc(Epetra_SerialDenseMatrix & matrix_X, Epetra_SerialDenseMatrix & xg, unsigned int & gp){
        Epetra_SerialDenseVector t(3), normal(3);
        Epetra_SerialDenseMatrix sigma(3,3), d_shape_functions(Mesh->face_type,2), dxi_matrix_x(3,2);
        sigma = manufacturedStress(xg(0,gp),xg(1,gp),xg(2,gp));
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            d_shape_functions(inode,0) = Mesh->D1_N_tri(gp,inode);
            d_shape_functions(inode,1) = Mesh->D2_N_tri(gp,inode);
        }
        dxi_matrix_x.Multiply('N','N',1.0,matrix_X,d_shape_functions,0.0);
        normal(0) = dxi_matrix_x(1,0)*dxi_matrix_x(2,1) - dxi_matrix_x(2,0)*dxi_matrix_x(1,1);
        normal(1) = dxi_matrix_x(2,0)*dxi_matrix_x(0,1) - dxi_matrix_x(0,0)*dxi_matrix_x(2,1);
        normal(2) = dxi_matrix_x(0,0)*dxi_matrix_x(1,1) - dxi_matrix_x(1,0)*dxi_matrix_x(0,1);
        if (xg(2,gp)!=0.0){
            normal.Scale(-1.0);
        }
        t.Multiply('N','N',1.0,sigma,normal,0.0);
        return t;
    }
    Epetra_SerialDenseVector get_forcing(double & x1, double & x2, double & x3, unsigned int & e_lid, unsigned int & gp){
        Epetra_SerialDenseVector f(3);
        Epetra_SerialDenseMatrix C(6,6);
        get_elasticity_tensor(e_lid,gp,C);
        double k = std::sqrt(2.0);
        double C11 = C(0,0); double C12 = C(0,1); double C13 = C(0,2); double C14 = C(0,3); double C15 = C(0,4); double C16 = C(0,5);
                             double C22 = C(1,1); double C23 = C(1,2); double C24 = C(1,3); double C25 = C(1,4); double C26 = C(1,5);
                                                  double C33 = C(2,2); double C34 = C(2,3); double C35 = C(2,4); double C36 = C(2,5);
                                                                       double C44 = C(3,3); double C45 = C(3,4); double C46 = C(3,5);
                                                                                            double C55 = C(4,4); double C56 = C(4,5);
                                                                                                                 double C66 = C(5,5);
        f(0) = 12*C11*x1*x1 + 6*k*C26*x2*x2 + (3*C56)/5 + (k*C14)/5 + k*C15;
        f(1) = 6*k*C16*x1*x1 + 12*C22*x2*x2 + (3*C46)/5 + C56 + (k*C25)/5;
        f(2) = 6*k*C15*x1*x1 + 6*k*C24*x2*x2 + (2*C45)/5 + C55 + (2*k*C36)/5;
        f.Scale(-1.0);
        return f;
    }
    
    void solve(bool doprint){
        Epetra_FECrsMatrix linearOperator(Copy,*FEGraph);
        Epetra_FEVector    rhs(*StandardMap);
        Epetra_Vector      lhs(*StandardMap);
        rhs.PutScalar(0.0);
        double dummy = 0.0;
        assembleMixedDirichletNeumann_inhomogeneousForcing(linearOperator,rhs);
        apply_dirichlet_conditions(linearOperator,rhs,dummy);
        aztecSolver(linearOperator,rhs,lhs,*Krylov);
        if (doprint){
            print_solution(lhs,"/Users/brian/Documents/GitHub/Trilinos_results/cee530/manufactured/manufactured.mtx");
        }
    }
    
    void setup_dirichlet_conditions(){
        n_bc_dof = 0;
        double x;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            x = Mesh->nodes_coord[3*node+1];
            if(x==0.0 || x==25.0/1000.0){
                n_bc_dof+=3;
            }
        }
        
        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x = Mesh->nodes_coord[3*node+1];
            if(x==0.0 || x==25.0/1000.0){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
        }
    }
    
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        Epetra_MultiVector v(*StandardMap,true);
        v.PutScalar(0.0);
        
        int node;
        double x, y, z;
        Epetra_SerialDenseVector u(3);
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            if(y==0.0 || y==25.0/1000.0){
                u = manufacturedSolution(x,y,z);
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
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            if(y==0.0 || y==25.0/1000.0){
                F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
                F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
                F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }
    
    void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix){
        int e_gid = Mesh->local_cells[e_lid];
        int n_gauss_cells = Mesh->n_gauss_cells;
        double c1 = 144.8969;
        double c2 = 14.2500;
        double c3 = 5.8442;
        double c4 = 7.5462;
        double c5 = 12.5580;
        transverse_isotropic_matrix(tangent_matrix,c1,c2,c3,c4,c5);
    }
    
    void transverse_isotropic_matrix(Epetra_SerialDenseMatrix & C, double & c1, double & c2, double & c3, double & c4, double & c5){
        C(0,0) = c1/16.0 + (9.0*c2)/32.0 + (9.0*c4)/32.0 + (3.0*c5)/8.0 + (3.0*std::sqrt(2.0)*c3)/16.0;
        C(0,1) = (3.0*c1)/16.0 + (3.0*c2)/32.0 + (3.0*c4)/32.0 - (3.0*c5)/8.0 + (5.0*std::sqrt(2.0)*c3)/16.0;
        C(0,2) = (3.0*c2)/8.0 - (3.0*c4)/8.0 + (std::sqrt(2.0)*c3)/8.0;
        C(0,3) = 0.0;
        C(0,4) = 0.0;
        C(0,5) = std::sqrt(6)*(c1/16.0 - (3.0*c2)/32.0 - (3.0*c4)/32.0 + c5/8.0 + (std::sqrt(2.0)*c3)/16.0);

        C(1,0) = (3.0*c1)/16.0 + (3.0*c2)/32.0 + (3.0*c4)/32.0 - (3.0*c5)/8.0 + (5*std::sqrt(2.0)*c3)/16.0;
        C(1,1) = (9.0*c1)/16.0 + c2/32.0 + c4/32.0 + (3.0*c5)/8.0 + (3.0*std::sqrt(2.0)*c3)/16.0;
        C(1,2) = c2/8.0 - c4/8.0 + (3.0*std::sqrt(2.0)*c3)/8.0;
        C(1,3) = 0.0;
        C(1,4) = 0.0;
        C(1,5) = -std::sqrt(6.0)*(c2/32.0 - (3.0*c1)/16.0 + c4/32.0 + c5/8.0 + (std::sqrt(2.0)*c3)/16.0);
        
        C(2,0) = (3.0*c2)/8.0 - (3.0*c4)/8.0 + std::sqrt(2.0)*c3/8.0;
        C(2,1) = c2/8.0 - c4/8.0 + (3.0*std::sqrt(2.0)*c3)/8.0;
        C(2,2) = c2/2.0 + c4/2.0;
        C(2,3) = 0.0;
        C(2,4) = 0.0;
        C(2,5) = (std::sqrt(3.0)*(2.0*c3 - std::sqrt(2.0)*c2 + std::sqrt(2.0)*c4))/8.0;
        
        C(3,0) = 0.0;
        C(3,1) = 0.0;
        C(3,2) = 0.0;
        C(3,3) = c4/4.0 + (3.0*c5)/4.0;
        C(3,4) = -(std::sqrt(3.0)*(c4 - c5))/4.0;
        C(3,5) = 0.0;
        
        C(4,0) = 0.0;
        C(4,1) = 0.0;
        C(4,2) = 0.0;
        C(4,3) = -(std::sqrt(3.0)*(c4 - c5))/4.0;
        C(4,4) = (3.0*c4)/4.0 + c5/4.0;
        C(4,5) = 0.0;
        
        C(5,0) = std::sqrt(6.0)*(c1/16.0 - (3.0*c2)/32.0 - (3.0*c4)/32.0 + c5/8.0 + std::sqrt(2.0)*c3/16.0);
        C(5,1) = -std::sqrt(6.0)*(c2/32.0 - (3.0*c1)/16.0 + c4/32.0 + c5/8.0 + (std::sqrt(2.0)*c3)/16.0);
        C(5,2) = (std::sqrt(3.0)*(2.0*c3 - std::sqrt(2.0)*c2 + std::sqrt(2.0)*c4))/8.0;
        C(5,3) = 0.0;
        C(5,4) = 0.0;
        C(5,5) = (3.0*c1)/8.0 + (3.0*c2)/16.0 + (3.0*c4)/16.0 + c5/4.0 - (3.0*std::sqrt(2.0)*c3)/8.0;
    }
    
    void get_elasticity_tensor_for_recovery(unsigned int & e_lid, Epetra_SerialDenseMatrix & tangent_matrix){
    }
};


#endif
