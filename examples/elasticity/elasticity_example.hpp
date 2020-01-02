/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef LINEARPATCHTEST_HPP
#define LINEARPATCHTEST_HPP

#include "linearizedElasticity.hpp"

class elasticity_example : public linearizedElasticity
{
public:

    Teuchos::ParameterList * Krylov;

    double E; // = 210000.0;
    double nu; // = 0.3;
    double lambda; // = E*nu/((1.0+nu)*(1.0-2.0*nu));
    double mu; // = E/(2.0*(1.0+nu));
    double kappa; // = lambda + 2.0*mu/3.0;
    double dummy;

    elasticity_example(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){

        Krylov = &Parameters.sublist("Krylov");

        E  = Teuchos::getParameter<double>(Parameters.sublist("Behavior"), "young");
        nu = Teuchos::getParameter<double>(Parameters.sublist("Behavior"), "poisson");
        dummy = Teuchos::getParameter<double>(Parameters.sublist("Behavior"), "factor");

        lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
        mu = E/(2.0*(1.0+nu));
        kappa = lambda + 2.0*mu/3.0;

        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh = new mesh(comm, mesh_file, 1.0);
        Comm = Mesh->Comm;

        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();

        setup_dirichlet_conditions();
    }

    ~elasticity_example(){
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

    void solve(bool doprint){
        Epetra_FECrsMatrix linearOperator(Copy,*FEGraph);
        Epetra_FEVector    rhs(*StandardMap);
        Epetra_Vector      lhs(*StandardMap);
        rhs.PutScalar(0.0);;
        assemblePureDirichlet_homogeneousForcing(linearOperator);
        apply_dirichlet_conditions(linearOperator,rhs,dummy);
        aztecSolver(linearOperator,rhs,lhs,*Krylov);
        if (doprint){
            print_solution(lhs,"/Users/brian/Documents/GitHub/TrilinosUQComp/results/elasticity_example/plate2020.mtx");
        }
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
            if(y==0.0||y==10.0){
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
            if(y==0.0||y==10.0){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
        }
    }

    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & factor){
        Epetra_MultiVector v(*StandardMap,true);
        v.PutScalar(0.0);

        int node;
        double x,y,z;
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            if (y==0.0||y==10.0){
                v[0][StandardMap->LID(3*node+0)] = 0.0;
                v[0][StandardMap->LID(3*node+1)] = factor*y;
                v[0][StandardMap->LID(3*node+2)] = 0.0;
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
            if (y==0.0||y==10.0){
                F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
                F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
                F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }

    void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tgm){
        double c1 = lambda+2.0*mu;
        tgm(0,0) = c1;     tgm(0,1) = lambda; tgm(0,2) = lambda; tgm(0,3) = 0.0;     tgm(0,4) = 0.0;    tgm(0,5) = 0.0;
        tgm(1,0) = lambda; tgm(1,1) = c1;     tgm(1,2) = lambda; tgm(1,3) = 0.0;     tgm(1,4) = 0.0;    tgm(1,5) = 0.0;
        tgm(2,0) = lambda; tgm(2,1) = lambda; tgm(2,2) = c1;     tgm(2,3) = 0.0;     tgm(2,4) = 0.0;    tgm(2,5) = 0.0;
        tgm(3,0) = 0.0;    tgm(3,1) = 0.0;    tgm(3,2) = 0.0;    tgm(3,3) = 2.0*mu;  tgm(3,4) = 0.0;    tgm(3,5) = 0.0;
        tgm(4,0) = 0.0;    tgm(4,1) = 0.0;    tgm(4,2) = 0.0;    tgm(4,3) = 0.0;     tgm(4,4) = 2.0*mu; tgm(4,5) = 0.0;
        tgm(5,0) = 0.0;    tgm(5,1) = 0.0;    tgm(5,2) = 0.0;    tgm(5,3) = 0.0;     tgm(5,4) = 0.0;    tgm(5,5) = 2.0*mu;
    }

    void get_elasticity_tensor_for_recovery(unsigned int & e_lid, Epetra_SerialDenseMatrix & tangent_matrix){
    }
};


#endif
