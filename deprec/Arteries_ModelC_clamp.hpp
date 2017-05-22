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
        
        for (unsigned int e=0; e<Mesh->n_cells/18; ++e){
            for (unsigned int j=0; j<12; ++j){
                phase.push_back(0);
            }
            for (unsigned int k=0; k<6; ++k){
                phase.push_back(1);
            }
        }
        
    }
    
    ~Interface_arteries(){
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
        //if (phase[e_gid]==0){
           /* mu1 = 6.8703/1000.0;
            mu2 = 0.71095/1000.0;
            mu4 = 95.051/1000.0;
            beta4 = 1401.2;
            theta = 0.80606;
            mu3 = 8.325/1000.0;
            beta3 = 2.2575;*/

            mu1 = 4.1543/1000.0;
            mu2 = 2.5084/1000.0;
            mu4 = 19.285/1000.0;
            beta4 = 500.02;
            theta = 0.80764;
            mu3 = 9.7227/1000.0;
            beta3 = 3.6537;
        //}
	/*if (phase[e_gid]==1){
	    mu1 = 2.3846/1000.0;
            mu2 = 2.478/1000.0;
            mu4 = 0.0026725/1000.0;
            beta4 = 48.069;
            theta = 0.63156;
            mu3 = 9.4331/1000.0;
            beta3 = 3.3953;
	}*/
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
        
        double I1 = C(0) + C(1) + C(2);
        double II1 = C(0)*C(0) + C(1)*C(1) + C(2)*C(2) + 2.0*C(3)*C(3) + 2.0*C(4)*C(4) + 2.0*C(5)*C(5);
        double I2 = (1.0/2.0)*(I1*I1-II1);
        double I4_1 = C(0)*M1(0) + C(1)*M1(1) + C(2)*M1(2) + 2.0*C(5)*M1(5) + 2.0*C(4)*M1(4) + 2.0*C(3)*M1(3);
        double I4_2 = C(0)*M2(0) + C(1)*M2(1) + C(2)*M2(2) + 2.0*C(5)*M2(5) + 2.0*C(4)*M2(4) + 2.0*C(3)*M2(3);
        double pI2 = std::sqrt(I2);
        
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
            scalarAB = (8.0*mu4 + 16.0*mu4*beta4*S4_1)*exp(beta4*S4_1);
            tensor_product(scalarAB,M1,M1,tangent_piola_isc,1.0);
        }
        if (I4_2>1.0){
            piola_isc += piola_ani2;
            scalarAB = (8.0*mu4 + 16.0*mu4*beta4*S4_2)*exp(beta4*S4_2);
            tensor_product(scalarAB,M2,M2,tangent_piola_isc,1.0);
        }
        
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
    
    double mu1;
    double mu2;
    double mu3;
    double mu4;
    double beta3;
    double beta4;
    double theta;
    
    std::vector<int> phase;
    
    Epetra_SerialDenseVector a;
    Epetra_SerialDenseVector b;
    
};

#endif
