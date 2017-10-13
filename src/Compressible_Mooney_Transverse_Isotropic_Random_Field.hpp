#ifndef COMPRESSIBLE_MOONEY_TRANSVERSE_ISOTROPIC_RANDOM_FIELD_HPP
#define COMPRESSIBLE_MOONEY_TRANSVERSE_ISOTROPIC_RANDOM_FIELD_HPP

#include "tensor_calculus.hpp"
#include "shinozukapp_layeredcomp_2d.hpp"
#include "hyperelasticity_setup_pp.hpp"

class TIMooney_RandomField : public hyperelasticity_setup
{
public:
    
    Teuchos::RCP<shinozuka_layeredcomp_2d> GRF_Generator;
    
    Epetra_SerialDenseVector mu1rf, mu2rf, mu3rf, mu4rf, mu5rf;
    double mu1, mu2, mu3, mu4, mu5;
    double mean_mu1, mean_mu2, mean_mu3, mean_mu4, mean_mu5, mu, trm;
    double beta3, beta4, beta5, ptrmbeta4, ptrmbeta5;
    double plyagl, cos_plyagl, sin_plyagl;
    int n_ply;
    std::vector<int> phase;
    
    TIMooney_RandomField(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
        
        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        n_ply = Teuchos::getParameter<int>(Parameters.sublist("Mesh"), "n_ply");
        Mesh = new mesh(comm, mesh_file);
        Comm = Mesh->Comm;
        
        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();
        
        setup_dirichlet_conditions();
        for (unsigned int e=0; e<Mesh->n_cells/32; ++e){
            for (unsigned int j=0; j<32; ++j){
                phase.push_back(j);
            }
        }
        
        double plyagldeg = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"angle");
        mean_mu1 = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu1");
        mean_mu2 = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu2");
        mean_mu3 = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu3");
        mean_mu4 = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu4");
        mean_mu5 = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"mu5");
        beta3    = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta3");
        beta4    = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta4");
        beta5    = Teuchos::getParameter<double>(paramList->sublist("TIMooney"),"beta5");
        plyagl   = 2.0*M_PI*plyagldeg/360.0;
        for (unsigned int i=0; i<5; i++){
            parameters(i) = 1.0e3*parameters(i);
        }
        int order = Teuchos::getParameter<int>(Parameters.sublist("Shinozuka"), "order");
        double L1 = Teuchos::getParameter<double>(Parameters.sublist("Shinozuka"), "lx");
        double L2 = Teuchos::getParameter<double>(Parameters.sublist("Shinozuka"), "ly");
        
        GRF_Generator = Teuchos::rcp(new shinozuka_layeredcomp_2d(order,L1,L2));
    }
    
    ~TIMooney_RandomField(){
    }
    
    void RandomFieldGenerator(int * seeds){
        
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
        
        GRF_Generator->rng.seed(seeds[0]);
        GRF_Generator->generator_gauss_points(w1_shino,*Mesh,phase);
        
        GRF_Generator->rng.seed(seeds[1]);
        GRF_Generator->generator_gauss_points(w2_shino,*Mesh,phase);
        
        GRF_Generator->rng.seed(seeds[2]);
        GRF_Generator->generator_gauss_points(w3_shino,*Mesh,phase);
        
        GRF_Generator->rng.seed(seeds[3]);
        GRF_Generator->generator_gauss_points(w4_shino,*Mesh,phase);
        
        GRF_Generator->rng.seed(seeds[4]);
        GRF_Generator->generator_gauss_points(w5_shino,*Mesh,phase);
        
        double delta = 0.10;
        double alpha = 1.0/(delta*delta);
        double beta;
        for (unsigned int i=0; i<w1_shino.Length(); ++i){
            beta = mean_mu1*delta*delta; mu1rf(i) = icdf_gamma(w1_shino(i),alpha,beta);
            beta = mean_mu2*delta*delta; mu2rf(i) = icdf_gamma(w2_shino(i),alpha,beta);
            beta = mean_mu3*delta*delta; mu3rf(i) = icdf_gamma(w3_shino(i),alpha,beta);
            beta = mean_mu4*delta*delta; mu4rf(i) = icdf_gamma(w4_shino(i),alpha,beta);
            beta = mean_mu5*delta*delta; mu5rf(i) = icdf_gamma(w5_shino(i),alpha,beta);
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
            if(coord==25.0){
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
            if (coord==25.0){
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
            if (coord==25.0){
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
            if (coord==25.0){
                F[0][StandardMap->LID(3*node+0)]   = 0.0;
                F[0][StandardMap->LID(3*node+dof)] = displacement;
                F[0][StandardMap->LID(3*node+2)]   = 0.0;
            }
        }
        //}
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
            sin_plyagl = -std::sin(plyagl);
        }
        else{
            cos_plyagl = std::cos(plyagl);
            sin_plyagl = -std::sin(plyagl);
        }
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
    
};

#endif