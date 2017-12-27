#ifndef COMPRESSIBLE_MOONEY_TRANSVERSE_ISOTROPIC_RANDOM_FIELD_HPP
#define COMPRESSIBLE_MOONEY_TRANSVERSE_ISOTROPIC_RANDOM_FIELD_HPP

#include "tensor_calculus.hpp"
#include "shinozukapp_layeredcomp_2d.hpp"
#include "compressibleHyperelasticity.hpp"

class TIMooney_RandomField : public compressibleHyperelasticity
{
public:

    Teuchos::RCP<shinozuka_layeredcomp_2d> GRF_Generator;

    Epetra_SerialDenseVector mu1rf, mu2rf, mu3rf, mu4rf, mu5rf;
    Epetra_SerialDenseVector omega;
    Epetra_SerialDenseVector mean_mu;
    double mu1, mu2, mu3, mu4, mu5;
    double mean_mu1, mean_mu2, mean_mu3, mean_mu4, mean_mu5, mu, trm;
    double beta3, beta4, beta5, ptrmbeta4, ptrmbeta5;
    double plyagl, cos_plyagl, sin_plyagl;
    double topcoord;
    std::vector<int> phase;

    TIMooney_RandomField(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){

        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh                  = new mesh(comm, mesh_file, 1.0);
        Comm                  = Mesh->Comm;

        StandardMap           = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap            = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap    = new Epetra_Import(*OverlapMap,*StandardMap);

        create_FECrsGraph();

        findtop();
        setup_dirichlet_conditions();
        for (unsigned int e=0; e<Mesh->n_cells/32; ++e){
            for (unsigned int j=0; j<32; ++j){
                phase.push_back(j);
            }
        }
        mean_mu.Resize(5);
        omega.Resize(6);
        int order = Teuchos::getParameter<int>(Parameters.sublist("Shinozuka"), "order");
        omega(4) = Teuchos::getParameter<double>(Parameters.sublist("Shinozuka"), "lx");
        omega(5) = Teuchos::getParameter<double>(Parameters.sublist("Shinozuka"), "ly");

        GRF_Generator = Teuchos::rcp(new shinozuka_layeredcomp_2d(order,omega(4),omega(5)));
    }

    ~TIMooney_RandomField(){
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

    void set_plyagl(double & Plyagl){
        plyagl = Plyagl;
    }

    void setParameters(Epetra_SerialDenseVector & parameters, Epetra_SerialDenseVector & exponents, Epetra_SerialDenseVector & hyperParameters){
        mean_mu = parameters;
        beta3   = -0.5;
        beta4   = exponents(0);
        beta5   = exponents(1);
        omega   = hyperParameters;
    }

    void RandomFieldGenerator(Epetra_IntSerialDenseVector & seeds){

        Epetra_SerialDenseVector w1_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);
        Epetra_SerialDenseVector w2_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);
        Epetra_SerialDenseVector w3_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);
        Epetra_SerialDenseVector w4_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);
        Epetra_SerialDenseVector w5_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);

        mu1rf.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);
        mu2rf.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);
        mu3rf.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);
        mu4rf.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);
        mu5rf.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);

        GRF_Generator->l1 = omega(4);
        GRF_Generator->l2 = omega(5);

        GRF_Generator->rng.seed(seeds(0));
        GRF_Generator->generator_gauss_points(w1_shino,*Mesh,phase);

        GRF_Generator->rng.seed(seeds(1));
        GRF_Generator->generator_gauss_points(w2_shino,*Mesh,phase);

        GRF_Generator->rng.seed(seeds(2));
        GRF_Generator->generator_gauss_points(w3_shino,*Mesh,phase);

        GRF_Generator->rng.seed(seeds(3));
        GRF_Generator->generator_gauss_points(w4_shino,*Mesh,phase);

        GRF_Generator->rng.seed(seeds(4));
        GRF_Generator->generator_gauss_points(w5_shino,*Mesh,phase);

        double alpha, beta;
        for (unsigned int i=0; i<w1_shino.Length(); ++i){
            alpha = 1.0/(omega(0)*omega(0)); beta = mean_mu(0)*omega(0)*omega(0);
            mu1rf(i) = icdf_gamma(w1_shino(i),alpha,beta);

            alpha = 1.0/(omega(1)*omega(1)); beta = mean_mu(1)*omega(1)*omega(1);
            mu2rf(i) = icdf_gamma(w2_shino(i),alpha,beta);

            alpha = 1.0/(omega(2)*omega(2)); beta = mean_mu(2)*omega(2)*omega(2);
            mu3rf(i) = icdf_gamma(w3_shino(i),alpha,beta);

            alpha = 1.0/(omega(3)*omega(3)); beta = mean_mu(3)*omega(3)*omega(3);
            mu4rf(i) = icdf_gamma(w4_shino(i),alpha,beta);

            alpha = (2.0*alpha)-1.0; beta = mean_mu(4)/alpha;
            mu5rf(i) = icdf_gamma(w5_shino(i),alpha,beta);
        }

    }

    double icdf_gamma(double & w, double & alpha, double & beta){
        double erfx = boost::math::erf<double>(w);
        double y = (1.0/2.0)*(1.0 + erfx);
        double yinv = boost::math::gamma_p_inv<double,double>(alpha,y);
        double z = yinv*beta;
        return z;
    }

    void get_matrix_and_rhs(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
        assemblePureDirichlet_homogeneousForcing(x,K,F);
    }

    void setup_dirichlet_conditions(){
        n_bc_dof = 0;
        double coord;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            coord = Mesh->nodes_coord[3*node+1];
            if(coord==0.0){
                n_bc_dof+=3;
            }
            if(coord==topcoord){
                n_bc_dof+=3;
            }
        }

        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+1];
            if (coord==0.0){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
            if (coord==topcoord){
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
        double coord;
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+1];
            if (coord==topcoord){
                v[0][StandardMap->LID(3*node+1)] = displacement;
            }
        }

        Epetra_MultiVector rhs_dir(*StandardMap,true);
        K.Apply(v,rhs_dir);
        F.Update(-1.0,rhs_dir,1.0);

        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+1];
            if (coord==0.0){
                F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
                F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
                F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
            }
            if (coord==topcoord){
                F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
                F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
                F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }

    void get_material_parameters(unsigned int & e_lid, unsigned int & gp){
        int e_gid = Mesh->local_cells[e_lid];
        int n_gauss_cells = Mesh->n_gauss_cells;
        mu1 = mu1rf(e_lid*n_gauss_cells+gp); mu2 = mu2rf(e_lid*n_gauss_cells+gp); mu3 = mu3rf(e_lid*n_gauss_cells+gp); mu4 = mu4rf(e_lid*n_gauss_cells+gp); mu5 = mu5rf(e_lid*n_gauss_cells+gp);
        mu = 2.0*mu1 + 4.0*mu2 + 2.0*mu3;
        trm = mu4 + 2.0*mu5;
        ptrmbeta4 = std::pow(trm,beta4);
        ptrmbeta5 = std::pow(trm,beta5);
        if (phase[e_gid] % 2){
            cos_plyagl = std::cos(plyagl);
            sin_plyagl = std::sin(plyagl);
        }
        else{
            cos_plyagl = std::cos(plyagl);
            sin_plyagl = -std::sin(plyagl);
        }
    }

    void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola){

        double det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);

        Epetra_SerialDenseMatrix C(3,3), CC(3,3), ddJ5(6,6);
        Epetra_SerialDenseVector LML(6), eye(6), dJ5(6), M(6), L(6), c(6);

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

        double I1  = C(0,0) + C(1,1) + C(2,2);
        double II1 = C(0,0)*C(0,0) + C(1,1)*C(1,1) + C(2,2)*C(2,2) + 2.0*C(1,2)*C(1,2) + 2.0*C(0,2)*C(0,2) + 2.0*C(0,1)*C(0,1);
        double I2  = (1.0/2.0)*(I1*I1-II1);
        double I3  = det*det;
        double I4  = C(0,0)*M(0) + C(1,1)*M(1) + C(2,2)*M(2) + 2.0*C(0,1)*M(5) + 2.0*C(0,2)*M(4) + 2.0*C(1,2)*M(3);
        double I5  = CC(0,0)*M(0) + CC(1,1)*M(1) + CC(2,2)*M(2) + 2.0*CC(0,1)*M(5) + 2.0*CC(0,2)*M(4) + 2.0*CC(1,2)*M(3);
        double J5  = I5 - I1*I4 + I2*trm;
        double pI3 = std::pow(I3,-beta3);
        double pI4 = std::pow(I4, beta4);
        double pJ5 = std::pow(J5, beta5);

        for (unsigned int i=0; i<6; ++i){
            dJ5(i) = J5*L(i) - I3*LML(i);
            piola_stress(i) = 2.0*mu1*eye(i) + 2.0*mu2*(I1*eye(i)-c(i)) + (2.0*mu3*det*det-mu)*L(i)
                            + (2.0/ptrmbeta4)*pI4*M(i) + (2.0/ptrmbeta5)*pJ5*dJ5(i) - 2.0*trm*pI3*L(i);
        }

        tensor_product(4.0*mu2,eye,eye,tangent_piola,0.0);
        sym_tensor_product(-4.0*mu2,eye,eye,tangent_piola,1.0);

        tensor_product(4.0*mu3*det*det,L,L,tangent_piola,1.0);
        sym_tensor_product(-4.0*mu3*det*det+2.0*mu,L,L,tangent_piola,1.0);

        tensor_product(J5,L,L,ddJ5,0.0);
        sym_tensor_product(-J5,L,L,ddJ5,1.0);

        tensor_product(-I3,L,LML,ddJ5,1.0);
        tensor_product(-I3,LML,L,ddJ5,1.0);

        sym_tensor_product(I3,L,LML,ddJ5,1.0);
        sym_tensor_product(I3,LML,L,ddJ5,1.0);

        tensor_product(4.0*beta4*pI4/(I4*ptrmbeta4),M,M,tangent_piola,1.0);
        tensor_product(4.0*beta5*pJ5/(J5*ptrmbeta5),dJ5,dJ5,tangent_piola,1.0);

        tensor_product(4.0*trm*beta3*pI3,L,L,tangent_piola,1.0);
        sym_tensor_product(4.0*trm*pI3,L,L,tangent_piola,1.0);

        ddJ5.Scale(4.0*pJ5/ptrmbeta5);
        tangent_piola += ddJ5;
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

    void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress){
    }

};

#endif
