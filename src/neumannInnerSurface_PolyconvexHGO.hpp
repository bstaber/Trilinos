/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef NEUMANNINNERSURFACE_POLYCONVEXHGO_HPP
#define NEUMANNINNERSURFACE_POLYCONVEXHGO_HPP

#include "tensor_calculus.hpp"
#include "laplacepp.hpp"
#include "nearlyIncompressibleHyperelasticity.hpp"

class neumannInnerSurface_PolyconvexHGO : public nearlyIncompressibleHyperelasticity
{
public:

    laplace * Laplace;

    double mu1, mu2, mu3, mu4, beta3, beta4, theta;
    Epetra_SerialDenseVector a, b;

    neumannInnerSurface_PolyconvexHGO(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){

        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        std::string boundary_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "boundary_file");
        unsigned int number_physical_groups = Teuchos::getParameter<unsigned int>(Parameters.sublist("Mesh"), "nb_phys_groups");
        std::string select_model = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "model");

        mu1   = Teuchos::getParameter<double>(Parameters.sublist(select_model), "mu1");
        mu2   = Teuchos::getParameter<double>(Parameters.sublist(select_model), "mu2");
        mu3   = Teuchos::getParameter<double>(Parameters.sublist(select_model), "mu3");
        mu4   = Teuchos::getParameter<double>(Parameters.sublist(select_model), "mu4");
        beta3 = Teuchos::getParameter<double>(Parameters.sublist(select_model), "beta3");
        beta4 = Teuchos::getParameter<double>(Parameters.sublist(select_model), "beta4");
        theta = Teuchos::getParameter<double>(Parameters.sublist(select_model), "theta");

        Mesh = new mesh(comm, Parameters); //mesh_file, 1000.0);
        Mesh->read_boundary_file(boundary_file,number_physical_groups);
        Comm = Mesh->Comm;

        Laplace = new laplace(*Mesh);

        //$$Laplace
        Epetra_FECrsMatrix matrix(Copy,*Laplace->FEGraph);
        Epetra_Vector psi(*Laplace->StandardMap);
        Epetra_Vector phi(*Laplace->StandardMap);
        Epetra_FEVector rhs(*Laplace->StandardMap);
        int bc_indx[2];
        double bc_val[2];
        bc_val[0] = 0.0; bc_val[1] = 1.0;
        //Problem number one with aztec
        bc_indx[0] = 2; bc_indx[1] = 3;
        Laplace->solve_aztec(Parameters.sublist("Laplace"), matrix, phi, rhs, &bc_indx[0], &bc_val[0]);
        Laplace->print_solution(phi, "laplace_inlet_to_outlet_aztec.mtx");
        //Problem number one with aztec
        bc_indx[0] = 0; bc_indx[1] = 1;
        Laplace->solve_aztec(Parameters.sublist("Laplace"), matrix, psi, rhs, &bc_indx[0], &bc_val[0]);
        Laplace->print_solution(psi, "laplace_inner_to_outer_aztec.mtx");
        //Get local directions
        Laplace->compute_local_directions(phi, psi);
        Laplace->compute_center_local_directions(phi, psi);
        //$$End

        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();

        a.Resize(3);
        b.Resize(3);
        setup_dirichlet_conditions();
    }

    ~neumannInnerSurface_PolyconvexHGO(){
    }

    void get_matrix_and_rhs(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
        assembleMixedDirichletDeformationDependentNeumann_homogeneousForcing(x,K,F);
    }

    void setup_dirichlet_conditions(){
        //setup_clamp();
        setup_semislipbc();
        //setup_slipbc_and_clamp();
        //setup_slipbc();
    }

    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        //apply_clamp(K,F,displacement);
        apply_semiselipbc(K,F,displacement);
        //apply_slipbc_and_clamp(K,F,displacement);
        //apply_slipbc(K,F,displacement);
    }

    void get_material_parameters(unsigned int & e_lid, unsigned int & gp){

        int n_gauss_points = Mesh->n_gauss_cells;
        int e_gid = Mesh->local_cells[e_lid];

        for (int i=0; i<3; ++i){
            a(i) = cos(theta)*Laplace->laplace_direction_one(n_gauss_points*e_lid+gp,i)
                 + sin(theta)*Laplace->laplace_direction_two_cross_one(n_gauss_points*e_lid+gp,i);
            b(i) = cos(theta)*Laplace->laplace_direction_one(n_gauss_points*e_lid+gp,i)
                 - sin(theta)*Laplace->laplace_direction_two_cross_one(n_gauss_points*e_lid+gp,i);
        }

    }

    void get_constitutive_tensors_static_condensation(Epetra_SerialDenseMatrix & deformation_gradient,
                                                      double & det, Epetra_SerialDenseVector & inverse_cauchy,
                                                      Epetra_SerialDenseVector & piola_isc,
                                                      Epetra_SerialDenseVector & piola_vol,
                                                      Epetra_SerialDenseMatrix & tangent_piola_isc,
                                                      Epetra_SerialDenseMatrix & tangent_piola_vol){

        model_C(deformation_gradient, det, inverse_cauchy, piola_isc, piola_vol, tangent_piola_isc, tangent_piola_vol);

    }

    void get_internal_pressure(double & theta, double & pressure, double & dpressure){
        double ptheta = std::pow(theta,beta3);
        pressure  = mu3*beta3*( (ptheta/theta) - (1.0/(ptheta*theta)) );
        dpressure = mu3*beta3*( (beta3-1.0)*(ptheta/(theta*theta)) + (beta3+1.0)/(ptheta*theta*theta) );
    }

    void get_material_parameters_for_recover(unsigned int & e_lid){

      for (int i=0; i<3; ++i){
          a(i) = cos(theta)*Laplace->laplace_direction_one_center(e_lid,i)
               + sin(theta)*Laplace->laplace_direction_two_cross_one_center(e_lid,i);
          b(i) = cos(theta)*Laplace->laplace_direction_one_center(e_lid,i)
               - sin(theta)*Laplace->laplace_direction_two_cross_one_center(e_lid,i);
      }

    }

    void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress){

        det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)
            - deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)
            - deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)
            + deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)
            + deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)
            - deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);

        double alpha = std::pow(det,-2.0/3.0);
        double beta = 1.0/(det*det);
        Epetra_SerialDenseMatrix eye(3,3);
        Epetra_SerialDenseMatrix M1(3,3), M2(3,3);
        Epetra_SerialDenseMatrix C(3,3), L(3,3);
        Epetra_SerialDenseMatrix piola_ani1(3,3), piola_ani2(3,3);

        eye(0,0) = 1.0; eye(0,1) = 0.0; eye(0,2) = 0.0;
        eye(1,0) = 0.0; eye(1,1) = 1.0; eye(1,2) = 0.0;
        eye(2,0) = 0.0; eye(2,1) = 0.0; eye(2,2) = 1.0;

        M1.Multiply('N','T',1.0,a,a,0.0);
        M2.Multiply('N','T',1.0,b,b,0.0);

        C.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);

        L(0,0) = (1.0/(det*det))*(C(1,1)*C(2,2)-C(1,2)*C(2,1));
        L(1,1) = (1.0/(det*det))*(C(0,0)*C(2,2)-C(0,2)*C(2,0));
        L(2,2) = (1.0/(det*det))*(C(0,0)*C(1,1)-C(0,1)*C(1,0));
        L(1,2) = (1.0/(det*det))*(C(0,2)*C(1,0)-C(0,0)*C(1,2));
        L(0,2) = (1.0/(det*det))*(C(0,1)*C(1,2)-C(0,2)*C(1,1));
        L(0,1) = (1.0/(det*det))*(C(0,2)*C(2,1)-C(0,1)*C(2,2));
        L(2,1) = L(1,2); L(2,0) = L(0,2); L(1,0) = L(0,1);

        double I1   = C(0,0) + C(1,1) + C(2,2);
        double II1  = C(0,0)*C(0,0) + C(1,1)*C(1,1) + C(2,2)*C(2,2) + 2.0*C(1,2)*C(1,2) + 2.0*C(0,2)*C(0,2) + 2.0*C(0,1)*C(0,1);
        double I2   = (1.0/2.0)*(I1*I1-II1);
        double I4_1 = C(0,0)*M1(0,0) + C(1,1)*M1(1,1) + C(2,2)*M1(2,2) + 2.0*C(0,1)*M1(0,1) + 2.0*C(0,2)*M1(0,2) + 2.0*C(1,2)*M1(1,2);
        double I4_2 = C(0,0)*M2(0,0) + C(1,1)*M2(1,1) + C(2,2)*M2(2,2) + 2.0*C(0,1)*M2(0,1) + 2.0*C(0,2)*M2(0,2) + 2.0*C(1,2)*M2(1,2);
        double pI2  = std::sqrt(I2);

        double S4_1 = (I4_1-1.0)*(I4_1-1.0);
        double S4_2 = (I4_2-1.0)*(I4_2-1.0);

        double ptheta   = std::pow(det,beta3);
        double pressure = mu3*beta3*( (ptheta/det) - (1.0/(ptheta*det)) );

        for (unsigned int i=0; i<3; ++i){
            for (unsigned int j=0; j<3; ++j){
                piola_stress(i,j) = 2.0*mu1*alpha*(eye(i,j)-(1.0/3.0)*L(i,j))
                + mu2*beta*( 3.0*pI2*(I1*eye(i,j)-C(i,j)) - 2.0*I2*pI2*L(i,j) )
                + det*pressure*L(i,j);
                piola_ani1(i,j) = 4.0*mu4*(I4_1-1.0)*exp(beta4*S4_1)*M1(i,j);
                piola_ani2(i,j) = 4.0*mu4*(I4_2-1.0)*exp(beta4*S4_2)*M2(i,j);
            }
        }

        if (I4_1>1.0){
            piola_stress += piola_ani1;
        }
        if (I4_2>1.0){
            piola_stress += piola_ani2;
        }

    }

    void setup_clamp(){
        n_bc_dof = 0;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            if (Mesh->nodes_to_boundaries(i,2)==1 || Mesh->nodes_to_boundaries(i,3)==1){
                n_bc_dof+=3;
            }
        }

        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            if (Mesh->nodes_to_boundaries(inode,2)==1 || Mesh->nodes_to_boundaries(inode,3)==1){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
        }
    }

    void setup_slipbc(){
        const int nodesblk_gid0 = 480-1;
        const int nodesblk_gid1 = 538-1;
        const int nodesblk_gid2 = 577-1;

        int node;
        n_bc_dof = 0;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            if (Mesh->nodes_to_boundaries(i,2)==1 || Mesh->nodes_to_boundaries(i,3)==1){
                node = Mesh->local_nodes[i];
                switch (node){
                    case nodesblk_gid0:
                        n_bc_dof+=2;
                        break;
                    case nodesblk_gid1:
                        n_bc_dof+=2;
                        break;
                    case nodesblk_gid2:
                        n_bc_dof+=2;
                        break;
                    default:
                        n_bc_dof+=1;
                        break;
                };
            }
        }

        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            if (Mesh->nodes_to_boundaries(inode,2)==1 || Mesh->nodes_to_boundaries(inode,3)==1){
                node = Mesh->local_nodes[inode];
                switch (node){
                    case nodesblk_gid0:
                        dof_on_boundary[indbc+0] = 3*inode+0;
                        dof_on_boundary[indbc+1] = 3*inode+1;
                        indbc+=2;
                        break;
                    case nodesblk_gid1:
                        dof_on_boundary[indbc+0] = 3*inode+0;
                        dof_on_boundary[indbc+1] = 3*inode+2;
                        indbc+=2;
                        break;
                    case nodesblk_gid2:
                        dof_on_boundary[indbc+0] = 3*inode+0;
                        dof_on_boundary[indbc+1] = 3*inode+2;
                        indbc+=2;
                        break;
                    default:
                        dof_on_boundary[indbc+0] = 3*inode+0;
                        indbc+=1;
                        break;
                };
            }
        }
    }

    void setup_semislipbc(){
        n_bc_dof = 0;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            if (Mesh->nodes_to_boundaries(i,2)==1){
                n_bc_dof+=1;
            }
            if (Mesh->nodes_to_boundaries(i,3)==1){
                n_bc_dof+=3;
            }
        }

        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            if (Mesh->nodes_to_boundaries(inode,2)==1){
                dof_on_boundary[indbc] = 3*inode+0;
                indbc+=1;
            }
            if (Mesh->nodes_to_boundaries(inode,3)==1){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
        }
    }

    void setup_slipbc_and_clamp(){
        Epetra_IntSerialDenseVector nodesblk_gid(3);
        nodesblk_gid(0) = 480-1;
        nodesblk_gid(1) = 481-1;
        nodesblk_gid(2) = 479-1;

        int node;
        n_bc_dof = 0;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            if (Mesh->nodes_to_boundaries(i,2)==1 || Mesh->nodes_to_boundaries(i,3)==1){
                node = Mesh->local_nodes[i];
                if (node==nodesblk_gid(0)||node==nodesblk_gid(1)||node==nodesblk_gid(2)){
                    n_bc_dof+=3;
                }
                else{
                    n_bc_dof+=1;
                }
            }
        }

        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            if (Mesh->nodes_to_boundaries(inode,2)==1 || Mesh->nodes_to_boundaries(inode,3)==1){
                node = Mesh->local_nodes[inode];
                if (node==nodesblk_gid(0)||node==nodesblk_gid(1)||node==nodesblk_gid(2)){
                    dof_on_boundary[indbc+0] = 3*inode+0;
                    dof_on_boundary[indbc+1] = 3*inode+1;
                    dof_on_boundary[indbc+2] = 3*inode+2;
                    indbc+=3;
                }
                else{
                    dof_on_boundary[indbc+0] = 3*inode+0;
                    indbc+=1;
                }
            }
        }
    }

    void apply_clamp(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        if (n_bc_dof>0){
            int node;
            for (int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
                node = Mesh->local_nodes[inode];
                if (Mesh->nodes_to_boundaries(inode,2)==1 || Mesh->nodes_to_boundaries(inode,3)==1){
                    F[0][StandardMap->LID(3*node+0)] = displacement;
                    F[0][StandardMap->LID(3*node+1)] = displacement;
                    F[0][StandardMap->LID(3*node+2)] = displacement;
                }
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }

    void apply_slipbc(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        const int nodesblk_gid0 = 480-1;
        const int nodesblk_gid1 = 538-1;
        const int nodesblk_gid2 = 577-1;

        if (n_bc_dof>0){
            int node;
            for (int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
                node = Mesh->local_nodes[inode];
                if (Mesh->nodes_to_boundaries(inode,2)==1 || Mesh->nodes_to_boundaries(inode,3)==1){
                    switch (node){
                        case nodesblk_gid0:
                            F[0][StandardMap->LID(3*node+0)] = displacement;
                            F[0][StandardMap->LID(3*node+1)] = displacement;
                            break;
                        case nodesblk_gid1:
                            F[0][StandardMap->LID(3*node+0)] = displacement;
                            F[0][StandardMap->LID(3*node+2)] = displacement;
                            break;
                        case nodesblk_gid2:
                            F[0][StandardMap->LID(3*node+0)] = displacement;
                            F[0][StandardMap->LID(3*node+2)] = displacement;
                            break;
                        default:
                            F[0][StandardMap->LID(3*node+0)] = displacement;
                            break;
                    };
                }
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }

    void apply_semiselipbc(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        if (n_bc_dof>0){
            int node;
            for (int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
                node = Mesh->local_nodes[inode];
                if (Mesh->nodes_to_boundaries(inode,2)==1){
                    F[0][StandardMap->LID(3*node+0)] = displacement;
                }
                if (Mesh->nodes_to_boundaries(inode,3)==1){
                    F[0][StandardMap->LID(3*node+0)] = displacement;
                    F[0][StandardMap->LID(3*node+1)] = displacement;
                    F[0][StandardMap->LID(3*node+2)] = displacement;
                }
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }

    void apply_slipbc_and_clamp(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        Epetra_IntSerialDenseVector nodesblk_gid(3);
        nodesblk_gid(0) = 480-1;
        nodesblk_gid(1) = 481-1;
        nodesblk_gid(2) = 479-1;

        if (n_bc_dof>0){
            int node;
            for (int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
                node = Mesh->local_nodes[inode];
                if (Mesh->nodes_to_boundaries(inode,2)==1 || Mesh->nodes_to_boundaries(inode,3)==1){
                    if (node==nodesblk_gid(0)||node==nodesblk_gid(1)||node==nodesblk_gid(2)){
                        F[0][StandardMap->LID(3*node+0)] = displacement;
                        F[0][StandardMap->LID(3*node+1)] = displacement;
                        F[0][StandardMap->LID(3*node+2)] = displacement;
                    }
                    else{
                        F[0][StandardMap->LID(3*node+0)] = displacement;
                    }
                }
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }

    void model_C(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseVector & L, Epetra_SerialDenseVector & piola_isc, Epetra_SerialDenseVector & piola_vol, Epetra_SerialDenseMatrix & tangent_piola_isc, Epetra_SerialDenseMatrix & tangent_piola_vol){

        det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)
            - deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)
            - deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)
            + deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)
            + deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)
            - deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);

        double alpha = std::pow(det,-2.0/3.0);

        Epetra_SerialDenseVector eye(6);
        Epetra_SerialDenseMatrix rightCauchy(3,3);
        Epetra_SerialDenseVector M1(6), M2(6);
        Epetra_SerialDenseVector C(6), D(6);
        Epetra_SerialDenseVector piola_nh(6), piola_ani1(6), piola_ani2(6);

        M1(0) = a(0)*a(0); M2(0) = b(0)*b(0);
        M1(1) = a(1)*a(1); M2(1) = b(1)*b(1);
        M1(2) = a(2)*a(2); M2(2) = b(2)*b(2);
        M1(3) = a(1)*a(2); M2(3) = b(1)*b(2);
        M1(4) = a(0)*a(2); M2(4) = b(0)*b(2);
        M1(5) = a(0)*a(1); M2(5) = b(0)*b(1);

        rightCauchy.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);

        eye(0) = 1.0; eye(1) = 1.0; eye(2) = 1.0; eye(3) = 0.0; eye(4) = 0.0; eye(5) = 0.0;

        C(0) = rightCauchy(0,0); C(1) = rightCauchy(1,1); C(2) = rightCauchy(2,2);
        C(3) = rightCauchy(1,2); C(4) = rightCauchy(0,2); C(5) = rightCauchy(0,1);

        L(0) = (1.0/(det*det))*(rightCauchy(1,1)*rightCauchy(2,2)-rightCauchy(1,2)*rightCauchy(2,1));
        L(1) = (1.0/(det*det))*(rightCauchy(0,0)*rightCauchy(2,2)-rightCauchy(0,2)*rightCauchy(2,0));
        L(2) = (1.0/(det*det))*(rightCauchy(0,0)*rightCauchy(1,1)-rightCauchy(0,1)*rightCauchy(1,0));
        L(3) = (1.0/(det*det))*(rightCauchy(0,2)*rightCauchy(1,0)-rightCauchy(0,0)*rightCauchy(1,2));
        L(4) = (1.0/(det*det))*(rightCauchy(0,1)*rightCauchy(1,2)-rightCauchy(0,2)*rightCauchy(1,1));
        L(5) = (1.0/(det*det))*(rightCauchy(0,2)*rightCauchy(2,1)-rightCauchy(0,1)*rightCauchy(2,2));

        double I1   = C(0) + C(1) + C(2);
        double II1  = C(0)*C(0) + C(1)*C(1) + C(2)*C(2) + 2.0*C(3)*C(3) + 2.0*C(4)*C(4) + 2.0*C(5)*C(5);
        double I2   = (1.0/2.0)*(I1*I1-II1);
        double I4_1 = C(0)*M1(0) + C(1)*M1(1) + C(2)*M1(2) + 2.0*C(5)*M1(5) + 2.0*C(4)*M1(4) + 2.0*C(3)*M1(3);
        double I4_2 = C(0)*M2(0) + C(1)*M2(1) + C(2)*M2(2) + 2.0*C(5)*M2(5) + 2.0*C(4)*M2(4) + 2.0*C(3)*M2(3);
        double pI2  = std::sqrt(I2);

        double S4_1 = (I4_1-1.0)*(I4_1-1.0);
        double S4_2 = (I4_2-1.0)*(I4_2-1.0);

        for (unsigned int i=0; i<6; ++i){
            D(i) = I1*eye(i) - C(i);
            piola_nh(i) = 2.0*alpha*mu1*(eye(i)-(1.0/3.0)*I1*L(i));
            piola_isc(i) = piola_nh(i) + (mu2/(det*det))*( 3.0*pI2*D(i) - 2.0*pI2*I2*L(i) );

            piola_ani1(i) = 4.0*mu4*(I4_1-1.0)*exp(beta4*S4_1)*M1(i);
            piola_ani2(i) = 4.0*mu4*(I4_2-1.0)*exp(beta4*S4_2)*M2(i);

            piola_vol(i) = det*L(i);
        }

        double scalarAB;

        //scalarAB = det;
        tensor_product(det,L,L,tangent_piola_vol,0.0);
        //scalarAB = -2.0*det;
        sym_tensor_product(-2.0*det,L,L,tangent_piola_vol,1.0);

        //scalarAB = -2.0/3.0;
        tensor_product(-2.0/3.0,piola_nh,L,tangent_piola_isc,0.0);
        tensor_product(-2.0/3.0,L,piola_nh,tangent_piola_isc,1.0);

        //scalarAB = -6.0*mu2*pI2/(det*det);
        tensor_product(-6.0*mu2*pI2/(det*det),D,eye,tangent_piola_isc,1.0);
        tensor_product(-6.0*mu2*pI2/(det*det),eye,D,tangent_piola_isc,1.0);

        //scalarAB = (-4.0/9.0)*mu1*alpha*I1 + 4.0*mu2*pI2*I2/(det*det);
        tensor_product((-4.0/9.0)*mu1*alpha*I1 + 4.0*mu2*pI2*I2/(det*det),L,L,tangent_piola_isc,1.0);

        //scalarAB = (4.0/3.0)*mu1*alpha*I1 + 4.0*mu2*pI2*I2/(det*det);
        sym_tensor_product((4.0/3.0)*mu1*alpha*I1 + 4.0*mu2*pI2*I2/(det*det),L,L,tangent_piola_isc,1.0);

        //scalarAB = 3.0*mu2/(det*det*pI2);
        tensor_product(3.0*mu2/(det*det*pI2),D,D,tangent_piola_isc,1.0);

        //scalarAB = 6.0*mu2*pI2/(det*det);
        tensor_product(6.0*mu2*pI2/(det*det),eye,eye,tangent_piola_isc,1.0);

        //scalarAB = -scalarAB;
        sym_tensor_product(-6.0*mu2*pI2/(det*det),eye,eye,tangent_piola_isc,1.0);

        if (I4_1>1.0){
            piola_isc += piola_ani1;
            //scalarAB = (8.0*mu4 + 16.0*mu4*beta4*S4_1)*exp(beta4*S4_1);
            tensor_product((8.0*mu4 + 16.0*mu4*beta4*S4_1)*exp(beta4*S4_1),M1,M1,tangent_piola_isc,1.0);
        }
        if (I4_2>1.0){
            piola_isc += piola_ani2;
            //scalarAB = (8.0*mu4 + 16.0*mu4*beta4*S4_2)*exp(beta4*S4_2);
            tensor_product((8.0*mu4 + 16.0*mu4*beta4*S4_2)*exp(beta4*S4_2),M2,M2,tangent_piola_isc,1.0);
        }
    }

    void model_B(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseVector & piola_isc, Epetra_SerialDenseVector & piola_vol, Epetra_SerialDenseMatrix & tangent_piola_isc, Epetra_SerialDenseMatrix & tangent_piola_vol){

        det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);

        double alpha = std::pow(det,-2.0/3.0);

        Epetra_SerialDenseMatrix rightCauchy(3,3);
        Epetra_SerialDenseVector M1(6), M2(6);
        Epetra_SerialDenseVector eye(6);
        Epetra_SerialDenseVector C(6), D(6), L(6);
        Epetra_SerialDenseVector dI4_1(6), dI4_2(6);
        Epetra_SerialDenseVector piola_nh(6), piola_ani1(6), piola_ani2(6);

        M1(0) = a(0)*a(0); M2(0) = b(0)*b(0);
        M1(1) = a(1)*a(1); M2(1) = b(1)*b(1);
        M1(2) = a(2)*a(2); M2(2) = b(2)*b(2);
        M1(3) = a(1)*a(2); M2(3) = b(1)*b(2);
        M1(4) = a(0)*a(2); M2(4) = b(0)*b(2);
        M1(5) = a(0)*a(1); M2(5) = b(0)*b(1);

        rightCauchy.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);

        eye(0) = 1.0; eye(1) = 1.0; eye(2) = 1.0; eye(3) = 0.0; eye(4) = 0.0; eye(5) = 0.0;

        C(0) = rightCauchy(0,0); C(1) = rightCauchy(1,1); C(2) = rightCauchy(2,2);
        C(3) = rightCauchy(1,2); C(4) = rightCauchy(0,2); C(5) = rightCauchy(0,1);

        L(0) = (1.0/(det*det))*(rightCauchy(1,1)*rightCauchy(2,2)-rightCauchy(1,2)*rightCauchy(2,1));
        L(1) = (1.0/(det*det))*(rightCauchy(0,0)*rightCauchy(2,2)-rightCauchy(0,2)*rightCauchy(2,0));
        L(2) = (1.0/(det*det))*(rightCauchy(0,0)*rightCauchy(1,1)-rightCauchy(0,1)*rightCauchy(1,0));
        L(3) = (1.0/(det*det))*(rightCauchy(0,2)*rightCauchy(1,0)-rightCauchy(0,0)*rightCauchy(1,2));
        L(4) = (1.0/(det*det))*(rightCauchy(0,1)*rightCauchy(1,2)-rightCauchy(0,2)*rightCauchy(1,1));
        L(5) = (1.0/(det*det))*(rightCauchy(0,2)*rightCauchy(2,1)-rightCauchy(0,1)*rightCauchy(2,2));

        double I1 = C(0) + C(1) + C(2);
        double II1 = C(0)*C(0) + C(1)*C(1) + C(2)*C(2) + 2.0*C(3)*C(3) + 2.0*C(4)*C(4) + 2.0*C(5)*C(5);
        double I2 = (1.0/2.0)*(I1*I1-II1);
        double I4_1 = C(0)*M1(0) + C(1)*M1(1) + C(2)*M1(2) + 2.0*C(5)*M1(5) + 2.0*C(4)*M1(4) + 2.0*C(3)*M1(3);
        double I4_2 = C(0)*M2(0) + C(1)*M2(1) + C(2)*M2(2) + 2.0*C(5)*M2(5) + 2.0*C(4)*M2(4) + 2.0*C(3)*M2(3);
        double pI2 = std::sqrt(I2);

        double S4_1 = (alpha*I4_1-1.0)*(alpha*I4_1-1.0);
        double S4_2 = (alpha*I4_2-1.0)*(alpha*I4_2-1.0);

        for (unsigned int i=0; i<6; ++i){
            D(i) = I1*eye(i) - C(i);
            dI4_1(i) = M1(i)-(1.0/3.0)*I4_1*L(i);
            dI4_2(i) = M2(i)-(1.0/3.0)*I4_2*L(i);
            piola_nh(i) = 2.0*alpha*mu1*(eye(i)-(1.0/3.0)*I1*L(i));
            piola_isc(i) = piola_nh(i) + (mu2/(det*det))*( 3.0*pI2*D(i) - 2.0*pI2*I2*L(i) );

            piola_ani1(i) = 4.0*mu4*alpha*(alpha*I4_1-1.0)*exp(beta4*S4_1)*(M1(i)-(1.0/3.0)*I4_1*L(i));
            piola_ani2(i) = 4.0*mu4*alpha*(alpha*I4_2-1.0)*exp(beta4*S4_2)*(M2(i)-(1.0/3.0)*I4_2*L(i));;

            piola_vol(i) = det*L(i);
        }

        double scalarAB;

        scalarAB = det;
        tensor_product(det,L,L,tangent_piola_vol,0.0);
        scalarAB = -2.0*det;
        sym_tensor_product(scalarAB,L,L,tangent_piola_vol,1.0);

        scalarAB = -2.0/3.0;
        tensor_product(scalarAB,piola_nh,L,tangent_piola_isc,0.0);
        tensor_product(scalarAB,L,piola_nh,tangent_piola_isc,1.0);

        scalarAB = -6.0*mu2*pI2/(det*det);
        tensor_product(scalarAB,D,eye,tangent_piola_isc,1.0);
        tensor_product(scalarAB,eye,D,tangent_piola_isc,1.0);

        scalarAB = (-4.0/9.0)*mu1*alpha*I1 + 4.0*mu2*pI2*I2/(det*det);
        tensor_product(scalarAB,L,L,tangent_piola_isc,1.0);

        scalarAB = (4.0/3.0)*mu1*alpha*I1 + 4.0*mu2*pI2*I2/(det*det);
        sym_tensor_product(scalarAB,L,L,tangent_piola_isc,1.0);

        scalarAB = 3.0*mu2/(det*det*pI2);
        tensor_product(scalarAB,D,D,tangent_piola_isc,1.0);

        scalarAB = 6.0*mu2*pI2/(det*det);
        tensor_product(scalarAB,eye,eye,tangent_piola_isc,1.0);

        scalarAB = -scalarAB;
        sym_tensor_product(scalarAB,eye,eye,tangent_piola_isc,1.0);

        if (I4_1>1.0){
            piola_isc += piola_ani1;
            scalarAB = -(8.0/3.0)*mu4*alpha*(alpha*I4_1-1.0)*exp(beta4*S4_1);
            tensor_product(scalarAB,dI4_1,L,tangent_piola_isc,1.0);
            scalarAB = (8.0/3.0)*mu4*alpha*alpha*exp(beta4*S4_1);
            tensor_product(scalarAB,dI4_1,dI4_1,tangent_piola_isc,1.0);
            scalarAB = 16.0*mu4*beta4*alpha*alpha*S4_1*exp(beta4*S4_1);
            tensor_product(scalarAB,dI4_1,dI4_1,tangent_piola_isc,1.0);
        }
        if (I4_2>1.0){
            piola_isc += piola_ani2;
            scalarAB = -(8.0/3.0)*mu4*alpha*(alpha*I4_2-1.0)*exp(beta4*S4_2);
            tensor_product(scalarAB,dI4_2,L,tangent_piola_isc,1.0);
            scalarAB = (8.0/3.0)*mu4*alpha*alpha*exp(beta4*S4_2);
            tensor_product(scalarAB,dI4_2,dI4_2,tangent_piola_isc,1.0);
            scalarAB = 16.0*mu4*beta4*alpha*alpha*S4_2*exp(beta4*S4_2);
            tensor_product(scalarAB,dI4_2,dI4_2,tangent_piola_isc,1.0);
        }

    }

    void model_A(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseVector & piola_isc, Epetra_SerialDenseVector & piola_vol, Epetra_SerialDenseMatrix & tangent_piola_isc, Epetra_SerialDenseMatrix & tangent_piola_vol){

        det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);

        double alpha = std::pow(det,-2.0/3.0);

        Epetra_SerialDenseMatrix rightCauchy(3,3);
        Epetra_SerialDenseVector M1(6), M2(6);
        Epetra_SerialDenseVector eye(6);
        Epetra_SerialDenseVector C(6), L(6), CC(6);
        Epetra_SerialDenseVector CMMC1(6), CMMC2(6);
        Epetra_SerialDenseVector D(6);
        Epetra_SerialDenseVector dK3_1(6), dK3_2(6);
        Epetra_SerialDenseVector piola_nh(6), piola_ani1(6), piola_ani2(6);

        M1(0) = a(0)*a(0);
        M1(1) = a(1)*a(1);
        M1(2) = a(2)*a(2);
        M1(3) = a(1)*a(2);
        M1(4) = a(0)*a(2);
        M1(5) = a(0)*a(1);

        M2(0) = b(0)*b(0);
        M2(1) = b(1)*b(1);
        M2(2) = b(2)*b(2);
        M2(3) = b(1)*b(2);
        M2(4) = b(0)*b(2);
        M2(5) = b(0)*b(1);

        rightCauchy.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);

        eye(0) = 1.0; eye(1) = 1.0; eye(2) = 1.0; eye(3) = 0.0; eye(4) = 0.0; eye(5) = 0.0;

        C(0) = rightCauchy(0,0); C(1) = rightCauchy(1,1); C(2) = rightCauchy(2,2);
        C(3) = rightCauchy(1,2); C(4) = rightCauchy(0,2); C(5) = rightCauchy(0,1);

        L(0) = (1.0/(det*det))*(rightCauchy(1,1)*rightCauchy(2,2)-rightCauchy(1,2)*rightCauchy(2,1));
        L(1) = (1.0/(det*det))*(rightCauchy(0,0)*rightCauchy(2,2)-rightCauchy(0,2)*rightCauchy(2,0));
        L(2) = (1.0/(det*det))*(rightCauchy(0,0)*rightCauchy(1,1)-rightCauchy(0,1)*rightCauchy(1,0));
        L(3) = (1.0/(det*det))*(rightCauchy(0,2)*rightCauchy(1,0)-rightCauchy(0,0)*rightCauchy(1,2));
        L(4) = (1.0/(det*det))*(rightCauchy(0,1)*rightCauchy(1,2)-rightCauchy(0,2)*rightCauchy(1,1));
        L(5) = (1.0/(det*det))*(rightCauchy(0,2)*rightCauchy(2,1)-rightCauchy(0,1)*rightCauchy(2,2));

        CMMC1(0) = 2.0*C(0)*M1(0) + 2.0*C(4)*M1(4) + 2.0*C(5)*M1(5);
        CMMC1(1) = 2.0*C(1)*M1(1) + 2.0*C(3)*M1(3) + 2.0*C(5)*M1(5);
        CMMC1(2) = 2.0*C(2)*M1(2) + 2.0*C(3)*M1(3) + 2.0*C(4)*M1(4);
        CMMC1(3) = C(1)*M1(3) + C(3)*M1(1) + C(2)*M1(3) + C(3)*M1(2) + C(4)*M1(5) + C(5)*M1(4);
        CMMC1(4) = C(0)*M1(4) + C(4)*M1(0) + C(2)*M1(4) + C(4)*M1(2) + C(3)*M1(5) + C(5)*M1(3);
        CMMC1(5) = C(0)*M1(5) + C(5)*M1(0) + C(1)*M1(5) + C(5)*M1(1) + C(3)*M1(4) + C(4)*M1(3);

        CMMC2(0) = 2.0*C(0)*M2(0) + 2.0*C(4)*M2(4) + 2.0*C(5)*M2(5);
        CMMC2(1) = 2.0*C(1)*M2(1) + 2.0*C(3)*M2(3) + 2.0*C(5)*M2(5);
        CMMC2(2) = 2.0*C(2)*M2(2) + 2.0*C(3)*M2(3) + 2.0*C(4)*M2(4);
        CMMC2(3) = C(1)*M2(3) + C(3)*M2(1) + C(2)*M2(3) + C(3)*M2(2) + C(4)*M2(5) + C(5)*M2(4);
        CMMC2(4) = C(0)*M2(4) + C(4)*M2(0) + C(2)*M2(4) + C(4)*M2(2) + C(3)*M2(5) + C(5)*M2(3);
        CMMC2(5) = C(0)*M2(5) + C(5)*M2(0) + C(1)*M2(5) + C(5)*M2(1) + C(3)*M2(4) + C(4)*M2(3);

        CC(0) = C(0)*C(0) + C(4)*C(4) + C(5)*C(5);
        CC(1) = C(1)*C(1) + C(3)*C(3) + C(5)*C(5);
        CC(2) = C(2)*C(2) + C(3)*C(3) + C(4)*C(4);
        CC(3) = C(1)*C(3) + C(2)*C(3) + C(4)*C(5);
        CC(4) = C(0)*C(4) + C(2)*C(4) + C(3)*C(5);
        CC(5) = C(0)*C(5) + C(1)*C(5) + C(3)*C(4);

        double I1 = C(0) + C(1) + C(2);
        double II1 = C(0)*C(0) + C(1)*C(1) + C(2)*C(2) + 2.0*C(3)*C(3) + 2.0*C(4)*C(4) + 2.0*C(5)*C(5);
        double I2 = (1.0/2.0)*(I1*I1-II1);
        double I4_1 = C(0)*M1(0) + C(1)*M1(1) + C(2)*M1(2) + 2.0*C(5)*M1(5) + 2.0*C(4)*M1(4) + 2.0*C(3)*M1(3);
        double I4_2 = C(0)*M2(0) + C(1)*M2(1) + C(2)*M2(2) + 2.0*C(5)*M2(5) + 2.0*C(4)*M2(4) + 2.0*C(3)*M2(3);
        double I5_1 = CC(0)*M1(0) + CC(1)*M1(1) + CC(2)*M1(2) + 2.0*CC(5)*M1(5) + 2.0*CC(4)*M1(4) + 2.0*CC(3)*M1(3);
        double I5_2 = CC(0)*M2(0) + CC(1)*M2(1) + CC(2)*M2(2) + 2.0*CC(5)*M2(5) + 2.0*CC(4)*M2(4) + 2.0*CC(3)*M2(3);
        double K3_1 = I1*I4_1-I5_1;
        double K3_2 = I1*I4_2-I5_2;
        double pI2 = std::sqrt(I2);

        for (unsigned int i=0; i<6; ++i){
            D(i) = I1*eye(i) - C(i);
            piola_nh(i) = 2.0*alpha*mu1*(eye(i)-(1.0/3.0)*I1*L(i));
            piola_isc(i) = piola_nh(i) + (mu2/(det*det))*( 3.0*pI2*D(i) - 2.0*pI2*I2*L(i) );

            dK3_1(i) = I4_1*eye(i) + I1*M1(i) - CMMC1(i);
            dK3_2(i) = I4_2*eye(i) + I1*M2(i) - CMMC2(i);
            piola_ani1(i) = 2.0*mu4*beta4*std::pow(K3_1-2.0,beta4-1.0)*dK3_1(i);
            piola_ani2(i) = 2.0*mu4*beta4*std::pow(K3_2-2.0,beta4-1.0)*dK3_2(i);

            piola_vol(i) = det*L(i);
        }

        double scalarAB;

        scalarAB = det;
        tensor_product(det,L,L,tangent_piola_vol,0.0);
        scalarAB = -2.0*det;
        sym_tensor_product(scalarAB,L,L,tangent_piola_vol,1.0);

        scalarAB = -2.0/3.0;
        tensor_product(scalarAB,piola_nh,L,tangent_piola_isc,0.0);
        tensor_product(scalarAB,L,piola_nh,tangent_piola_isc,1.0);

        scalarAB = -6.0*mu2*pI2/(det*det);
        tensor_product(scalarAB,D,eye,tangent_piola_isc,1.0);
        tensor_product(scalarAB,eye,D,tangent_piola_isc,1.0);

        scalarAB = (-4.0/9.0)*mu1*alpha*I1 + 4.0*mu2*pI2*I2/(det*det);
        tensor_product(scalarAB,L,L,tangent_piola_isc,1.0);

        scalarAB = (4.0/3.0)*mu1*alpha*I1 + 4.0*mu2*pI2*I2/(det*det);
        sym_tensor_product(scalarAB,L,L,tangent_piola_isc,1.0);

        scalarAB = 3.0*mu2/(det*det*pI2);
        tensor_product(scalarAB,D,D,tangent_piola_isc,1.0);

        scalarAB = 6.0*mu2*pI2/(det*det);
        tensor_product(scalarAB,eye,eye,tangent_piola_isc,1.0);

        scalarAB = -scalarAB;
        sym_tensor_product(scalarAB,eye,eye,tangent_piola_isc,1.0);

        if (K3_1>2.0){
            piola_isc += piola_ani1;
            scalarAB = 4.0*mu4*beta4*(beta4-1.0)*std::pow(K3_1-2.0,beta4-2.0);
            tensor_product(scalarAB,dK3_1,dK3_1,tangent_piola_isc,1.0);

            scalarAB = 4.0*mu4*beta4*std::pow(K3_1-2.0,beta4-1.0);
            tensor_product(scalarAB,eye,M1,tangent_piola_isc,1.0);
            tensor_product(scalarAB,M1,eye,tangent_piola_isc,1.0);
            scalarAB = -scalarAB;
            sym_tensor_product(scalarAB,M1,eye,tangent_piola_isc,1.0);
            sym_tensor_product(scalarAB,eye,M1,tangent_piola_isc,1.0);
        }
        if (K3_2>2.0){
            piola_isc += piola_ani2;
            scalarAB = 4.0*mu4*beta4*(beta4-1.0)*std::pow(K3_2-2.0,beta4-2.0);
            tensor_product(scalarAB,dK3_2,dK3_2,tangent_piola_isc,1.0);

            scalarAB = 4.0*mu4*beta4*std::pow(K3_2-2.0,beta4-1.0);
            tensor_product(scalarAB,eye,M2,tangent_piola_isc,1.0);
            tensor_product(scalarAB,M2,eye,tangent_piola_isc,1.0);
            scalarAB = -scalarAB;
            sym_tensor_product(scalarAB,M2,eye,tangent_piola_isc,1.0);
            sym_tensor_product(scalarAB,eye,M2,tangent_piola_isc,1.0);
        }
    }

};

#endif
