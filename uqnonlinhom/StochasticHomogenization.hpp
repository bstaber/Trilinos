#ifndef STOCHASTICHOMOGENIZATION_HPP
#define STOCHASTICHOMOGENIZATION_HPP

#include "tensor_calculus.hpp"
#include "compressibleHyperelasticity.hpp"

class StochasticHomogenization : public compressibleHyperelasticity
{
public:

    double p1m, p2m, p3m, sm;
    double p1f, p2f, p3f, sf;

    StochasticHomogenization(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){

        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh = new mesh(comm,mesh_file,1.0);
        Comm = Mesh->Comm;

        StandardMap        = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap         = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();

        setup_dirichlet_conditions();
    }

    ~StochasticHomogenization(){
    }

    void set_parameters(Epetra_SerialDenseVector & x){
        p1m = x(0);
        p2m = x(1);
        p3m = x(2);
        p1f = x(3);
        p2f = x(4);
        p3f = x(5);
        sm = 2.0*p1m * 4.0*p2m + 2.0*p3m;
        sf = 2.0*p1f * 4.0*p2f + 2.0*p3f;
    }

    void get_matrix_and_rhs(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
        assemblePureDirichlet_homogeneousForcing(x,K,F);
    }

    void setup_dirichlet_conditions(){
        n_bc_dof = 0;
        double x,y,z;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            x    = Mesh->nodes_coord[3*node+0];
            y    = Mesh->nodes_coord[3*node+1];
            x    = Mesh->nodes_coord[3*node+2];
            if(x==0.0 || y==0.0 || z==0.0 || x==1.0 || y==1.0 || z==1.0){
                n_bc_dof+=3;
            }
        }

        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x    = Mesh->nodes_coord[3*node+0];
            y    = Mesh->nodes_coord[3*node+1];
            x    = Mesh->nodes_coord[3*node+2];
            if(x==0.0 || y==0.0 || z==0.0 || x==1.0 || y==1.0 || z==1.0){
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
            x    = Mesh->nodes_coord[3*node+0];
            y    = Mesh->nodes_coord[3*node+1];
            x    = Mesh->nodes_coord[3*node+2];
            if(x==0.0 || y==0.0 || z==0.0 || x==1.0 || y==1.0 || z==1.0){
                v[0][StandardMap->LID(3*node+0)] = 0.0;
                v[0][StandardMap->LID(3*node+1)] = 0.0;
                v[0][StandardMap->LID(3*node+2)] = 0.0;
            }
        }

        Epetra_MultiVector rhs_dir(*StandardMap,true);
        K.Apply(v,rhs_dir);
        F.Update(-1.0,rhs_dir,1.0);

        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x    = Mesh->nodes_coord[3*node+0];
            y    = Mesh->nodes_coord[3*node+1];
            x    = Mesh->nodes_coord[3*node+2];
            if(x==0.0 || y==0.0 || z==0.0 || x==1.0 || y==1.0 || z==1.0){
                F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
                F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
                F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }

    void get_material_parameters(unsigned int & e_lid, unsigned int & gp){
      //todo
    }

    void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola){

        double det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);

        Epetra_SerialDenseMatrix C(3,3);
        Epetra_SerialDenseVector c(6);
        Epetra_SerialDenseVector eye(6), L(6);

        C.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);

        L(0) = (1.0/(det*det))*(C(1,1)*C(2,2)-C(1,2)*C(2,1));
        L(1) = (1.0/(det*det))*(C(0,0)*C(2,2)-C(0,2)*C(2,0));
        L(2) = (1.0/(det*det))*(C(0,0)*C(1,1)-C(0,1)*C(1,0));
        L(3) = (1.0/(det*det))*(C(0,2)*C(1,0)-C(0,0)*C(1,2));
        L(4) = (1.0/(det*det))*(C(0,1)*C(1,2)-C(0,2)*C(1,1));
        L(5) = (1.0/(det*det))*(C(0,2)*C(2,1)-C(0,1)*C(2,2));

        eye(0) = 1.0; eye(1) = 1.0; eye(2) = 1.0; eye(3) = 0.0; eye(4) = 0.0; eye(5) = 0.0;
        c(0) = C(0,0); c(1) = C(1,1); c(2) = C(2,2); c(3) = C(1,2); c(4) = C(0,2); c(5) = C(0,1);

        double I1         = C(0,0) + C(1,1) + C(2,2);
        double dpressure  = p3m*(det-1.0) - sm/det;

        for (unsigned int i=0; i<6; ++i){
            piola_stress(i) = 2.0*p1m*eye(i) + 2.0*p2m*(I1*eye(i)-c(i)) + det*dpressure*L(i);
        }

        double scalarAB = 4.0*p2m;
        tensor_product(scalarAB,eye,eye,tangent_piola,0.0);
        scalarAB = -4.0*p2m;
        sym_tensor_product(scalarAB,eye,eye,tangent_piola,1.0);
        scalarAB = p3m*det*(2.0*det-1.0);
        tensor_product(scalarAB,L,L,tangent_piola,1.0);
        scalarAB = -2.0*det*dpressure;
        sym_tensor_product(scalarAB,L,L,tangent_piola,1.0);
    }

    Epetra_SerialDenseVector get_neumannBc(Epetra_SerialDenseMatrix & matrix_X, Epetra_SerialDenseMatrix & xg, unsigned int & gp){
        Epetra_SerialDenseVector h(3);
        std::cout << "**Err: Not using that method in this example!\n";
        return h;
    }

    Epetra_SerialDenseVector get_forcing(double & x1, double & x2, double & x3, unsigned int & e_lid, unsigned int & gp){
        Epetra_SerialDenseVector f(3);
        std::cout << "**Err: Not using that method in this example!\n";
        return f;
    }

    void get_material_parameters_for_recover(unsigned int & e_lid){
        std::cout << "**Err: Not using that method in this example!\n";
    }

    void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress){
        std::cout << "**Err: Not using that method in this example!\n";
    }

};

#endif
