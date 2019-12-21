/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef LINEARPATCHTEST_HPP
#define LINEARPATCHTEST_HPP

#include "linearizedElasticity.hpp"

class linearPatchTest : public linearizedElasticity
{
public:

    Teuchos::ParameterList * Krylov;

    linearPatchTest(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){

        Krylov = &Parameters.sublist("Krylov");

        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh = new mesh(comm, mesh_file, 1000.0);
        Comm = Mesh->Comm;

        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();

        setup_dirichlet_conditions();
    }

    ~linearPatchTest(){
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

    double solve(bool doprint){
        Epetra_FECrsMatrix linearOperator(Copy,*FEGraph);
        Epetra_FEVector    rhs(*StandardMap);
        Epetra_Vector      lhs(*StandardMap);
        rhs.PutScalar(0.0);
        double dummy = 0.0;
        assemblePureDirichlet_homogeneousForcing(linearOperator);
        apply_dirichlet_conditions(linearOperator,rhs,dummy);
        aztecSolver(linearOperator,rhs,lhs,*Krylov);
        if (doprint){
            print_solution(lhs,"/Users/brian/Documents/GitHub/TrilinosUQComp/results/cee530/linearpatchtest/linearPatchTest2019.mtx");
        }
        double error = errorL2(lhs);
        return error;
    }

    double errorL2(Epetra_Vector & uStandardMap){

        Epetra_Vector u(*OverlapMap);
        u.Import(uStandardMap, *ImportToOverlapMap, Insert);
        double totalError;
        double error = 0.0;
        double normVH;
        double gauss_weight;
        int n_gauss_points = Mesh->n_gauss_cells;
        int e_gid, node;

        //Epetra_SerialDenseVector epsilon(6);
        //Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
        Epetra_SerialDenseVector uExact(3), vH(3);
        Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3), X_I(3,Mesh->el_type), u_I(3,Mesh->el_type);
        Epetra_SerialDenseMatrix u_G(3,n_gauss_points), x_G(3,n_gauss_points);

        for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
            e_gid = Mesh->local_cells[e_lid];
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
                X_I(0,inode) = Mesh->nodes_coord[3*node+0];
                X_I(1,inode) = Mesh->nodes_coord[3*node+1];
                X_I(2,inode) = Mesh->nodes_coord[3*node+2];
                u_I(0,inode) = u[OverlapMap->LID(3*node+0)];
                u_I(1,inode) = u[OverlapMap->LID(3*node+1)];
                u_I(2,inode) = u[OverlapMap->LID(3*node+2)];
            }
            x_G.Multiply('N','N',1.0,X_I,Mesh->N_cells,0.0);
            u_G.Multiply('N','N',1.0,u_I,Mesh->N_cells,0.0);
            for (unsigned int gp=0; gp<n_gauss_points; ++gp){
                gauss_weight = Mesh->gauss_weight_cells(gp);
                uExact = exactSolution(x_G(0,gp),x_G(1,gp),x_G(2,gp));
                vH(0) = uExact(0) - u_G(0,gp);
                vH(1) = uExact(1) - u_G(1,gp);
                vH(2) = uExact(2) - u_G(2,gp);
                normVH = vH.Norm2();
                /*for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                 dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,e_lid);
                 dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,e_lid);
                 dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,e_lid);
                 }
                 compute_B_matrices(dx_shape_functions,matrix_B);
                 epsilon.Multiply('N','N',1.0,matrix_B,vector_u,0.0);*/
                error += gauss_weight*normVH*normVH*Mesh->detJac_cells(e_lid,gp);
            }
        }
        Comm->SumAll(&error,&totalError,1);
        totalError = std::sqrt(totalError);
        return totalError;
    }

    Epetra_SerialDenseVector exactSolution(double & x1, double & x2, double & x3){
        Epetra_SerialDenseVector u(3);
        u(0) = x1; //0.1*x1+0.2*x2+0.4*x3;
        u(1) = 0.0; //0.4*x1+0.5*x2+0.1*x3;
        u(2) = 0.0; //0.05*x1+0.25*x2+0.65*x3;
        return u;
    }

    void setup_dirichlet_conditions(){
        n_bc_dof = 0;
        double x,y,z;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            if(x==0.0||y==0.0||z==0.0||x==1.0||y==1.0||z==1.0){
                n_bc_dof+=3;
            }
        }

        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            if(x==0.0||y==0.0||z==0.0||x==1.0||y==1.0||z==1.0){
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
        double x,y,z;
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            if (x==0.0||y==0.0||z==0.0||x==1.0||y==1.0||z==1.0){
                Epetra_SerialDenseVector u(3);
                u = exactSolution(x,y,z);
                v[0][StandardMap->LID(3*node+0)] = u(0); //0.1*(0.1*x + 0.08*y + 0.05*z + 0.04*x*y + 0.03*y*z + 0.08*z*x);
                v[0][StandardMap->LID(3*node+1)] = u(1); //x*y; //y*x; //0.1*(0.05*x + 0.04*y + 0.1*z + 0.07*x*y + 0.03*y*z + 0.08*z*x);
                v[0][StandardMap->LID(3*node+2)] = u(2); //x*y*z; //x*y*z; //0.1*(0.75*x + 0.09*y + 0.06*z + 0.07*x*y + 0.03*y*z + 0.08*z*x);
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
            if (x==0.0||y==0.0||z==0.0||x==1.0||y==1.0||z==1.0){
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
        double c1 = 144.8969*1.0e3;
        double c2 = 14.2500*1.0e3;
        double c3 = 5.8442*1.0e3;
        double c4 = 7.5462*1.0e3;
        double c5 = 12.5580*1.0e3;
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
