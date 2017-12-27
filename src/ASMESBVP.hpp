#ifndef ASMESBVP_HPP
#define ASMESBVP_HPP

#include "shinozukapp_layeredcomp.hpp"
#include "tensor_calculus.hpp"
#include "linearelasticity_setup_pp.hpp"

class ASMESBVP : public LinearizedElasticity
{
public:

    ASMESBVP(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){

        Krylov = &Parameters.sublist("Krylov");

        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        Mesh = new mesh(comm, mesh_file, 1.0);
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

        int order = Teuchos::getParameter<int>(Parameters.sublist("Shinozuka"), "order");
        double L1 = Teuchos::getParameter<double>(Parameters.sublist("Shinozuka"), "lx");
        double L2 = Teuchos::getParameter<double>(Parameters.sublist("Shinozuka"), "ly");
        double L3 = Teuchos::getParameter<double>(Parameters.sublist("Shinozuka"), "lz");

        GRF_Generator = Teuchos::rcp(new shinozuka_layeredcomp(order,L1,L2,L3));

        solution = new Epetra_Vector(*StandardMap);
    }

    ~ASMESBVP(){
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

    void solveOneRealization(double & bcDisp, int * seeds){

        Epetra_SerialDenseVector w1_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);
        Epetra_SerialDenseVector w2_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);
        Epetra_SerialDenseVector w3_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);
        Epetra_SerialDenseVector w4_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);
        Epetra_SerialDenseVector w5_shino(Mesh->n_local_cells*Mesh->n_gauss_cells);

        m1.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);
        m2.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);
        m3.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);
        m4.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);
        m5.Resize(Mesh->n_local_cells*Mesh->n_gauss_cells);

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

        double deltaN = 0.0970;
        double deltaM4 = 0.0461;
        double deltaM5 = 0.0952;
        double alpha; double beta = 1.0;
        double Psi1, Psi2;

        for (unsigned int i=0; i<w1_shino.Length(); ++i){

            alpha = 3.0/(2.0*deltaN*deltaN); beta = 1.0;
            Psi1 = icdf_gamma(w1_shino(i),alpha,beta);
            m1(i) = (deltaN*deltaN/3.0)*2.0*Psi1;

            alpha = 3.0/(2.0*deltaN*deltaN) - 1.0/2.0;
            Psi2 = icdf_gamma(w2_shino(i),alpha,beta);
            m2(i) = (deltaN*deltaN/3.0)*( 2.0*Psi2 + w3_shino(i)*w3_shino(i) );

            m3(i) = (deltaN*deltaN/3.0)*std::sqrt(2.0*Psi1)*w3_shino(i);

            alpha = 1.0/(deltaM4*deltaM4); beta = 1.0*deltaM4*deltaM4;
            m4(i) = icdf_gamma(w4_shino(i),alpha,beta);
            alpha = 1.0/(deltaM5*deltaM5); beta = 1.0*deltaM5*deltaM5;
            m5(i) = icdf_gamma(w5_shino(i),alpha,beta);

        }

        Epetra_FECrsMatrix linearOperator(Copy,*FEGraph);
        Epetra_FEVector    rhs(*StandardMap);
        Epetra_Vector      lhs(*StandardMap);

        rhs.PutScalar(0.0);
        lhs.PutScalar(0.0);

        assemblePureDirichlet_homogeneousForcing(linearOperator);
        apply_dirichlet_conditions(linearOperator,rhs,bcDisp);

        Epetra_LinearProblem problem;
        AztecOO solver;

        problem.SetOperator(&linearOperator);
        problem.SetLHS(&lhs);
        problem.SetRHS(&rhs);
        solver.SetProblem(problem);
        solver.SetParameters(*Krylov);

        solver.Iterate(2000,1e-6);

        *solution = lhs;

    }

    double icdf_gamma(double & w, double & alpha, double & beta){
        double erfx = boost::math::erf<double>(w);
        double y = (1.0/2.0)*(1.0 + erfx);
        double yinv = boost::math::gamma_p_inv<double,double>(alpha,y);
        double z = yinv*beta;
        return z;
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

    void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix){

        int e_gid = Mesh->local_cells[e_lid];
        int n_gauss_cells = Mesh->n_gauss_cells;

        double epsilon = 1.0e-6;
        Epetra_SerialDenseMatrix M(6,6);
        double M1 = m1(e_lid*n_gauss_cells+gp) + epsilon;
        double M2 = m2(e_lid*n_gauss_cells+gp) + epsilon;
        double M3 = m3(e_lid*n_gauss_cells+gp);
        double M4 = m4(e_lid*n_gauss_cells+gp) + epsilon;
        double M5 = m5(e_lid*n_gauss_cells+gp) + epsilon;

        transverse_isotropic_matrix(M,M1,M2,M3,M4,M5);

        double c1 = 144.8969*1.0e3;
        double c2 = 14.2500*1.0e3;
        double c3 = 5.8442*1.0e3;
        double c4 = 7.5462*1.0e3;
        double c5 = 12.5580*1.0e3;

        Epetra_SerialDenseMatrix sqrtmCmoy(6,6);
        double constant = std::sqrt(c1*c1-2.0*c1*c2+c2*c2+4.0*c3*c3);

        double d1 = (std::sqrt(2.0)*(c1-c2+constant)*std::sqrt(c1+c2+constant))/(4.0*constant) + (std::sqrt(2.0)*(c2-c1+constant)*std::sqrt(c1+c2-constant))/(4.0*constant);
        double d2 = (std::sqrt(2.0)*(c2-c1+constant)*std::sqrt(c1+c2+constant))/(4.0*constant) + (std::sqrt(2.0)*(c1-c2+constant)*std::sqrt(c1+c2-constant))/(4.0*constant);
        double d3 = std::sqrt(2.0)*c3*( std::sqrt(c1+c2+constant) - std::sqrt(c1+c2-constant) )/( 2.0*constant );
        double d4 = std::sqrt(c4);
        double d5 = std::sqrt(c5);

        transverse_isotropic_matrix(sqrtmCmoy,d1,d2,d3,d4,d5);

        Epetra_SerialDenseMatrix AtimesB(6,6);

        AtimesB.Multiply('N','N',1.0,M,sqrtmCmoy,0.0);
        tangent_matrix.Multiply('N','N',1.0,sqrtmCmoy,AtimesB,0.0);

        if(phase[e_gid] % 2){
            tangent_matrix(0,5) = -tangent_matrix(0,5);
            tangent_matrix(5,0) = -tangent_matrix(5,0);
            tangent_matrix(1,5) = -tangent_matrix(1,5);
            tangent_matrix(5,1) = -tangent_matrix(5,1);
            tangent_matrix(2,5) = -tangent_matrix(2,5);
            tangent_matrix(5,2) = -tangent_matrix(5,2);
            tangent_matrix(3,4) = -tangent_matrix(3,4);
            tangent_matrix(4,3) = -tangent_matrix(4,3);
        }

        tangent_matrix.Scale(1.0/(1.0+epsilon));

    }

    void recover_cauchy_stress(std::string & filename, int * seeds){

        Epetra_SerialDenseVector w1_shino(Mesh->n_local_cells);
        Epetra_SerialDenseVector w2_shino(Mesh->n_local_cells);
        Epetra_SerialDenseVector w3_shino(Mesh->n_local_cells);
        Epetra_SerialDenseVector w4_shino(Mesh->n_local_cells);
        Epetra_SerialDenseVector w5_shino(Mesh->n_local_cells);

        m1.Resize(Mesh->n_local_cells);
        m2.Resize(Mesh->n_local_cells);
        m3.Resize(Mesh->n_local_cells);
        m4.Resize(Mesh->n_local_cells);
        m5.Resize(Mesh->n_local_cells);

        double xi = 0.0; double eta = 0.0; double zeta = 0.0;

        GRF_Generator->rng.seed(seeds[0]);
        GRF_Generator->generator_one_gauss_point(w1_shino,*Mesh,phase,xi,eta,zeta);

        GRF_Generator->rng.seed(seeds[1]);
        GRF_Generator->generator_one_gauss_point(w2_shino,*Mesh,phase,xi,eta,zeta);

        GRF_Generator->rng.seed(seeds[2]);
        GRF_Generator->generator_one_gauss_point(w3_shino,*Mesh,phase,xi,eta,zeta);

        GRF_Generator->rng.seed(seeds[3]);
        GRF_Generator->generator_one_gauss_point(w4_shino,*Mesh,phase,xi,eta,zeta);

        GRF_Generator->rng.seed(seeds[4]);
        GRF_Generator->generator_one_gauss_point(w5_shino,*Mesh,phase,xi,eta,zeta);

        double deltaN = 0.0970;
        double deltaM4 = 0.0461;
        double deltaM5 = 0.0952;
        double alpha; double beta = 1.0;
        double Psi1, Psi2;

        for (unsigned int i=0; i<w1_shino.Length(); ++i){

            alpha = 3.0/(2.0*deltaN*deltaN); beta = 1.0;
            Psi1 = icdf_gamma(w1_shino(i),alpha,beta);
            m1(i) = (deltaN*deltaN/3.0)*2.0*Psi1;

            alpha = 3.0/(2.0*deltaN*deltaN) - 1.0/2.0;
            Psi2 = icdf_gamma(w2_shino(i),alpha,beta);
            m2(i) = (deltaN*deltaN/3.0)*( 2.0*Psi2 + w3_shino(i)*w3_shino(i) );

            m3(i) = (deltaN*deltaN/3.0)*std::sqrt(2.0*Psi1)*w3_shino(i);

            alpha = 1.0/(deltaM4*deltaM4); beta = 1.0*deltaM5*deltaM5;
            m4(i) = icdf_gamma(w4_shino(i),alpha,beta);
            alpha = 1.0/(deltaM4*deltaM4); beta = 1.0*deltaM5*deltaM5;
            m5(i) = icdf_gamma(w5_shino(i),alpha,beta);

        }
        compute_center_cauchy_stress(*solution, filename, true, true);
    }

    void get_elasticity_tensor_for_recovery(unsigned int & e_lid, Epetra_SerialDenseMatrix & tangent_matrix){

        int e_gid = Mesh->local_cells[e_lid];

        double epsilon = 1.0e-6;
        Epetra_SerialDenseMatrix M(6,6);
        double M1 = m1(e_lid) + epsilon;
        double M2 = m2(e_lid) + epsilon;
        double M3 = m3(e_lid);
        double M4 = m4(e_lid) + epsilon;
        double M5 = m5(e_lid) + epsilon;

        transverse_isotropic_matrix(M,M1,M2,M3,M4,M5);

        double c1 = 144.8969*1.0e3;
        double c2 = 14.2500*1.0e3;
        double c3 = 5.8442*1.0e3;
        double c4 = 7.5462*1.0e3;
        double c5 = 12.5580*1.0e3;

        Epetra_SerialDenseMatrix sqrtmCmoy(6,6);
        double constant = std::sqrt(c1*c1-2.0*c1*c2+c2*c2+4.0*c3*c3);

        double d1 = (std::sqrt(2.0)*(c1-c2+constant)*std::sqrt(c1+c2+constant))/(4.0*constant) + (std::sqrt(2.0)*(c2-c1+constant)*std::sqrt(c1+c2-constant))/(4.0*constant);
        double d2 = (std::sqrt(2.0)*(c2-c1+constant)*std::sqrt(c1+c2+constant))/(4.0*constant) + (std::sqrt(2.0)*(c1-c2+constant)*std::sqrt(c1+c2-constant))/(4.0*constant);
        double d3 = std::sqrt(2.0)*c3*( std::sqrt(c1+c2+constant) - std::sqrt(c1+c2-constant) )/( 2.0*constant );
        double d4 = std::sqrt(c4);
        double d5 = std::sqrt(c5);

        transverse_isotropic_matrix(sqrtmCmoy,d1,d2,d3,d4,d5);

        Epetra_SerialDenseMatrix AtimesB(6,6);

        AtimesB.Multiply('N','N',1.0,M,sqrtmCmoy,0.0);
        tangent_matrix.Multiply('N','N',1.0,sqrtmCmoy,AtimesB,0.0);

        if(phase[e_gid] % 2){
            tangent_matrix(0,5) = -tangent_matrix(0,5);
            tangent_matrix(5,0) = -tangent_matrix(5,0);
            tangent_matrix(1,5) = -tangent_matrix(1,5);
            tangent_matrix(5,1) = -tangent_matrix(5,1);
            tangent_matrix(2,5) = -tangent_matrix(2,5);
            tangent_matrix(5,2) = -tangent_matrix(5,2);
            tangent_matrix(3,4) = -tangent_matrix(3,4);
            tangent_matrix(4,3) = -tangent_matrix(4,3);
        }

        tangent_matrix.Scale(1.0/(1.0+epsilon));

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

    int print_solution(std::string fileName){

        int NumTargetElements = 0;
        if (Comm->MyPID()==0){
            NumTargetElements = 3*Mesh->n_nodes;
        }
        Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
        Epetra_Export ExportOnRoot(*StandardMap,MapOnRoot);
        Epetra_MultiVector lhs_root(MapOnRoot,true);
        lhs_root.Export(*solution,ExportOnRoot,Insert);

        int error = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(),lhs_root,0,0,false);

        return error;

    }

    Epetra_Vector * solution;
    Teuchos::ParameterList * Krylov;
    Teuchos::RCP<shinozuka_layeredcomp> GRF_Generator;

    Epetra_SerialDenseVector m1;
    Epetra_SerialDenseVector m2;
    Epetra_SerialDenseVector m3;
    Epetra_SerialDenseVector m4;
    Epetra_SerialDenseVector m5;

    std::vector<int> phase;

    double _deltaN;
    double _deltaM4;
    double _deltaM5;

};


#endif
