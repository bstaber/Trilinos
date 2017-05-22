#ifndef Interface_arteries_hpp
#define Interface_arteries_hpp

#include "tensor_calculus.hpp"
#include "hyperelasticity_setup_pp.hpp"

class Interface_arteries : public hyperelasticity_setup
{
public:
    
    Interface_arteries(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
        
        std::string mesh_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "mesh_file");
        std::string boundary_file = Teuchos::getParameter<std::string>(Parameters.sublist("Mesh"), "boundary_file");
        unsigned int number_physical_groups = Teuchos::getParameter<unsigned int>(Parameters.sublist("Mesh"), "nb_phys_groups");
        
        Mesh = new mesh(comm, mesh_file);
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
        //$$End
        
        StandardMap = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
        OverlapMap = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
        ImportToOverlapMap = new Epetra_Import(*OverlapMap,*StandardMap);
        create_FECrsGraph();
        
        a.Resize(3);
        b.Resize(3);
        setup_dirichlet_conditions();
        
        int count_med = 0;
        int count_adv = 0;
        for (unsigned int e=0; e<Mesh->n_cells/18; ++e){
            for (unsigned int j=0; j<12; ++j){
                phase.push_back(0);
                map_egid_elayer.push_back(count_med);
                count_med+=1;
            }
            for (unsigned int k=0; k<6; ++k){
                phase.push_back(1);
                map_egid_elayer.push_back(count_adv);
                count_adv+=1;
            }
        }
        
        N.Resize(4);
        xi.Resize(4); eta.Resize(4); zeta.Resize(4);
        double alpha = (5.0 - sqrt(5.0))/20.0;
        double beta = (5.0 + 3.0*sqrt(5.0))/20.0;
        xi[0] = alpha; eta[0] = alpha; zeta[0] = alpha;
        xi[1] = alpha; eta[1] = alpha; zeta[1] = beta;
        xi[2] = alpha; eta[2] = beta;  zeta[2] = alpha;
        xi[3] = beta;  eta[3] = alpha; zeta[3] = alpha;
    }
    
    ~Interface_arteries(){
    }
    
    void get_media_adventitia(unsigned int & n_cells_med, unsigned int & n_nodes_med, unsigned int & n_cells_adv, unsigned int & n_nodes_adv, std::string & path){
        
        std::ifstream connectivity_file_med, connectivity_file_adv;
        connectivity_file_med.open(path+"connectivity_med.txt");
        connectivity_file_adv.open(path+"connectivity_adv.txt");
        
        c1_med.Resize(n_nodes_med);
        c2_med.Resize(n_nodes_med);
        u1_med.Resize(n_nodes_med);
        mu4_med.Resize(n_nodes_med);
        c1_adv.Resize(n_nodes_adv);
        c2_adv.Resize(n_nodes_adv);
        u1_adv.Resize(n_nodes_adv);
        mu4_adv.Resize(n_nodes_adv);
        
        if (connectivity_file_med.is_open()){
            cells_nodes_med.Resize(4*n_cells_med);
            for (unsigned int e=0; e<4*n_cells_med; ++e){
                connectivity_file_med >> cells_nodes_med[e];
                cells_nodes_med[e] = cells_nodes_med[e]-1;
            }
            connectivity_file_med.close();
        }
        else{
            std::cout << "Couldn't open the connectivity file for the media.\n";
        }
        
        if (connectivity_file_adv.is_open()){
            cells_nodes_adv.Resize(4*n_cells_adv);
            for (unsigned int e=0; e<4*n_cells_adv; ++e){
                connectivity_file_adv >> cells_nodes_adv[e];
                cells_nodes_adv[e] = cells_nodes_adv[e]-1;
            }
            connectivity_file_adv.close();
        }
        else{
            std::cout << "Couldn't open the connectivity file for the adventitia.\n";
        }
    }
    
    void get_matrix_and_rhs(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
        assemble_dirichlet_neumann(x,K,F);
    }
    
    void setup_dirichlet_conditions(){
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
    
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
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
    
    void get_material_parameters(unsigned int & e_lid, unsigned int & gp){
        
        int n_gauss_points = Mesh->n_gauss_cells;
        int e_gid = Mesh->local_cells[e_lid];
        int e_layer = map_egid_elayer[e_gid];
        
        tetra4::shape_functions(N,xi[gp],eta[gp],zeta[gp]);
        
        if (phase[e_gid]==0){
            theta = 0.80521;
            beta3 = 4.7963;
            beta4 = 3000.0;
            
            c1 = N(0)*c1_med(cells_nodes_med(4*e_layer)) + N(1)*c1_med(cells_nodes_med(4*e_layer+1)) + N(2)*c1_med(cells_nodes_med(4*e_layer+2)) + N(4)*c1_med(cells_nodes_med(4*e_layer+3));
            c2 = N(0)*c2_med(cells_nodes_med(4*e_layer)) + N(1)*c2_med(cells_nodes_med(4*e_layer+1)) + N(2)*c2_med(cells_nodes_med(4*e_layer+2)) + N(4)*c2_med(cells_nodes_med(4*e_layer+3));
            u1 = N(0)*u1_med(cells_nodes_med(4*e_layer)) + N(1)*u1_med(cells_nodes_med(4*e_layer+1)) + N(2)*u1_med(cells_nodes_med(4*e_layer+2)) + N(4)*u1_med(cells_nodes_med(4*e_layer+3));
            mu4 = N(0)*mu4_med(cells_nodes_med(4*e_layer)) + N(1)*mu4_med(cells_nodes_med(4*e_layer+1)) + N(2)*mu4_med(cells_nodes_med(4*e_layer+2)) + N(4)*mu4_med(cells_nodes_med(4*e_layer+3));
            
            mu1 = (1.0/2.0)*c2*u1;
            mu2 = (1.0/(std::sqrt(3.0)*3.0))*c2*(1.0-u1);
            mu3 = c1/(2.0*beta3*beta3);
        }
        if (phase[e_gid]==1){
            theta = 0.61851;
            beta3 = 3.6252;
            beta4 = 51.062;
            
            c1 = N(0)*c1_adv(cells_nodes_adv(4*e_layer)) + N(1)*c1_adv(cells_nodes_adv(4*e_layer+1)) + N(2)*c1_adv(cells_nodes_adv(4*e_layer+2)) + N(3)*c1_adv(cells_nodes_adv(4*e_layer+3));
            c2 = N(0)*c2_adv(cells_nodes_adv(4*e_layer)) + N(1)*c2_adv(cells_nodes_adv(4*e_layer+1)) + N(2)*c2_adv(cells_nodes_adv(4*e_layer+2)) + N(4)*c2_adv(cells_nodes_adv(4*e_layer+3));
            u1 = N(0)*u1_adv(cells_nodes_adv(4*e_layer)) + N(1)*u1_adv(cells_nodes_adv(4*e_layer+1)) + N(2)*u1_adv(cells_nodes_adv(4*e_layer+2)) + N(4)*u1_adv(cells_nodes_adv(4*e_layer+3));
            mu4 = N(0)*mu4_adv(cells_nodes_adv(4*e_layer)) + N(1)*mu4_adv(cells_nodes_adv(4*e_layer+1)) + N(2)*mu4_adv(cells_nodes_adv(4*e_layer+2)) + N(4)*mu4_adv(cells_nodes_adv(4*e_layer+3));
            
            mu1 = (1.0/2.0)*c2*u1;
            mu2 = (1.0/(std::sqrt(3.0)*3.0))*c2*(1.0-u1);
            mu3 = c1/(2.0*beta3*beta3);
        }
        
        for (int i=0; i<3; ++i){
            a(i) = cos(theta)*Laplace->laplace_direction_one(n_gauss_points*e_lid+gp,i) + sin(theta)*Laplace->laplace_direction_two_cross_one(n_gauss_points*e_lid+gp,i);
            b(i) = cos(theta)*Laplace->laplace_direction_one(n_gauss_points*e_lid+gp,i) - sin(theta)*Laplace->laplace_direction_two_cross_one(n_gauss_points*e_lid+gp,i);
        }
        
    }
    
    void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector & piola_stress, Epetra_SerialDenseMatrix & tangent_piola){
        std::cerr << "**Err: Using static condensation method!\n";
    }
    
    void get_constitutive_tensors_static_condensation(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseVector & piola_isc, Epetra_SerialDenseVector & piola_vol, Epetra_SerialDenseMatrix & tangent_piola_isc, Epetra_SerialDenseMatrix & tangent_piola_vol){
        
        det = deformation_gradient(0,0)*deformation_gradient(1,1)*deformation_gradient(2,2)-deformation_gradient(0,0)*deformation_gradient(1,2)*deformation_gradient(2,1)-deformation_gradient(0,1)*deformation_gradient(1,0)*deformation_gradient(2,2)+deformation_gradient(0,1)*deformation_gradient(1,2)*deformation_gradient(2,0)+deformation_gradient(0,2)*deformation_gradient(1,0)*deformation_gradient(2,1)-deformation_gradient(0,2)*deformation_gradient(1,1)*deformation_gradient(2,0);
        
        double alpha = std::pow(det,-2.0/3.0);
        
        Epetra_SerialDenseMatrix rightCauchy(3,3);
        Epetra_SerialDenseVector M1(6), M2(6);
        Epetra_SerialDenseVector eye(6);
        Epetra_SerialDenseVector C(6);
        Epetra_SerialDenseVector D(6);
        Epetra_SerialDenseVector L(6);
        Epetra_SerialDenseVector piola_nh(6);
        
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
        
        double I1 = C(0) + C(1) + C(2);
        double II1 = C(0)*C(0) + C(1)*C(1) + C(2)*C(2) + 2.0*C(3)*C(3) + 2.0*C(4)*C(4) + 2.0*C(5)*C(5);
        double I2 = (1.0/2.0)*(I1*I1-II1);
        double I4_1 = C(0)*M1(0) + C(1)*M1(1) + C(2)*M1(2) + 2.0*C(5)*M1(5) + 2.0*C(4)*M1(4) + 2.0*C(3)*M1(3);
        double I4_2 = C(0)*M2(0) + C(1)*M2(1) + C(2)*M2(2) + 2.0*C(5)*M2(5) + 2.0*C(4)*M2(4) + 2.0*C(3)*M2(3);
        double pI2 = std::sqrt(I2);
        
        double S4_1 = 0.0;
        double T4_1 = 0.0;
        double U4_1 = 0.0;
        if (I4_1>1.0){
            U4_1 = 1.0;
            T4_1 = I4_1-1.0;
            S4_1 = T4_1*T4_1;
        }
        double S4_2 = 0.0;
        double T4_2 = 0.0;
        double U4_2 = 0.0;
        if (I4_2>1.0){
            U4_2 = 1.0;
            T4_2 = I4_2-1.0;
            S4_2 = T4_2*T4_2;
        }
        
        for (unsigned int i=0; i<6; ++i){
            D(i) = I1*eye(i) - C(i);
            piola_nh(i) = 2.0*alpha*mu1*(eye(i)-(1.0/3.0)*I1*L(i));
            piola_isc(i) = piola_nh(i) + (mu2/(det*det))*( 3.0*pI2*D(i) - 2.0*pI2*I2*L(i) ) + 4.0*mu4*T4_1*exp(beta4*S4_1)*M1(i) + 4.0*mu4*T4_2*exp(beta4*S4_2)*M2(i);
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
        
        scalarAB = (8.0*mu4*U4_1 + 16.0*mu4*beta4*S4_1)*exp(beta4*S4_1);
        tensor_product(scalarAB,M1,M1,tangent_piola_isc,1.0);
        
        scalarAB = (8.0*mu4*U4_2 + 16.0*mu4*beta4*S4_2)*exp(beta4*S4_2);
        tensor_product(scalarAB,M2,M2,tangent_piola_isc,1.0);
        
    }
    
    void get_internal_pressure(double & theta, double & pressure, double & dpressure){
        double ptheta = std::pow(theta,beta3);
        pressure = mu3*beta3*( (ptheta/theta) - (1.0/(ptheta*theta)) );
        dpressure = mu3*beta3*( (beta3-1.0)*(ptheta/(theta*theta)) + (beta3+1.0)/(ptheta*theta*theta) );
    }
    
    void get_material_parameters_for_recover(unsigned int & e_lid, double & xi, double & eta, double & zeta){
    }
    void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress){
    }
    
    laplace * Laplace;
    
    double c1;
    double c2;
    double u1;
    double mu1;
    double mu2;
    double mu3;
    double mu4;
    double beta3;
    double beta4;
    double theta;
    
    std::vector<int> phase;
    std::vector<int> map_egid_elayer;
    Epetra_IntSerialDenseVector cells_nodes_med;
    Epetra_IntSerialDenseVector cells_nodes_adv;
    Epetra_SerialDenseVector c1_med;
    Epetra_SerialDenseVector c2_med;
    Epetra_SerialDenseVector u1_med;
    Epetra_SerialDenseVector mu4_med;
    Epetra_SerialDenseVector c1_adv;
    Epetra_SerialDenseVector c2_adv;
    Epetra_SerialDenseVector u1_adv;
    Epetra_SerialDenseVector mu4_adv;
    
    Epetra_SerialDenseVector a;
    Epetra_SerialDenseVector b;
    Epetra_SerialDenseVector N;
    Epetra_SerialDenseVector xi;
    Epetra_SerialDenseVector eta;
    Epetra_SerialDenseVector zeta;
    
};

#endif
