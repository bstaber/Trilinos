#include "hyperelasticity_setup_pp.hpp"
#include "fepp.hpp"

hyperelasticity_setup::hyperelasticity_setup(){
}

hyperelasticity_setup::~hyperelasticity_setup(){
}

void hyperelasticity_setup::create_FECrsGraph(){
    
    FEGraph = new Epetra_FECrsGraph(Copy,*StandardMap,100);
    int eglob, node;
    int *index;
    index = new int [3*Mesh->el_type];
    
    for (int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        eglob = Mesh->local_cells[e_lid];
        for (int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
            for (int ddl=0; ddl<3; ++ddl){
                index[3*inode+ddl] = 3*node+ddl;
            }
        }
        for (int i=0; i<3*Mesh->el_type; ++i){
            for (int j=0; j<3*Mesh->el_type; ++j){
                FEGraph->InsertGlobalIndices(1, &index[i], 1, &index[j]);
            }
        }
        
    }
    Comm->Barrier();
    FEGraph->GlobalAssemble();
    delete[] index;
}

void hyperelasticity_setup::assemble_dirichlet(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){

    int error;
    
    F.PutScalar(0.0);
    K.PutScalar(0.0);
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    material_stiffness_and_rhs_dirichlet(u,K,F);
    
    Comm->Barrier();
    
    error=K.GlobalAssemble();
    error=K.FillComplete();
    error=F.GlobalAssemble();
}

void hyperelasticity_setup::assemble_dirichlet_static_condensation(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    
    F.PutScalar(0.0);
    K.PutScalar(0.0);
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    material_stiffness_and_rhs_static_condensation(u,K,F);
    
    Comm->Barrier();
    
    K.GlobalAssemble();
    K.FillComplete();
    F.GlobalAssemble();
}

void hyperelasticity_setup::assemble_dirichlet_live_neumann_static_condensation(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){

    F.PutScalar(0.0);
    K.PutScalar(0.0);
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    material_stiffness_and_rhs_static_condensation(u,K,F);
    force_stiffness_rhs_live_pressure(u,K,F);
    
    Comm->Barrier();
    K.GlobalAssemble();
    K.FillComplete();
    F.GlobalAssemble();
}

void hyperelasticity_setup::assemble_dirichlet_dead_neumann_static_condensation(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    
    F.PutScalar(0.0);
    K.PutScalar(0.0);
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    material_stiffness_and_rhs_static_condensation(u,K,F);
    force_dead_pressure(F);
    
    Comm->Barrier();
    K.GlobalAssemble();
    K.FillComplete();
    F.GlobalAssemble();
}

void hyperelasticity_setup::assemble_dirichlet_dead_neumann(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    
    F.PutScalar(0.0);
    K.PutScalar(0.0);
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    material_stiffness_and_rhs_dirichlet(u,K,F);
    force_dead_pressure(F);
    
    Comm->Barrier();
    
    K.GlobalAssemble();
    K.FillComplete();
    F.GlobalAssemble();
}

void hyperelasticity_setup::material_stiffness_and_rhs_dirichlet(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F){

    int node, e_gid, error;
    int n_gauss_points = Mesh->n_gauss_cells;
    double gauss_weight;

    int *Indices_tetra;
    Indices_tetra = new int [3*Mesh->el_type];
    
    Epetra_SerialDenseMatrix Ke(3*Mesh->el_type,3*Mesh->el_type);
    Epetra_SerialDenseVector Re(3*Mesh->el_type);
    
    Epetra_SerialDenseMatrix deformation_gradient(3,3);
    Epetra_SerialDenseMatrix tangent_matrix(6,6);
    Epetra_SerialDenseVector scd_piola_stress(6);
    Epetra_SerialDenseMatrix block_scd_piola_stress(9,9);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
    
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
    Epetra_SerialDenseMatrix B_times_TM(3*Mesh->el_type,6);
    Epetra_SerialDenseMatrix matrix_BG(9,3*Mesh->el_type);
    Epetra_SerialDenseMatrix BG_times_BSPS(3*Mesh->el_type,9);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            matrix_x(0,inode) = u[OverlapMap->LID(3*node+0)] + Mesh->nodes_coord[3*node+0];
            matrix_x(1,inode) = u[OverlapMap->LID(3*node+1)] + Mesh->nodes_coord[3*node+1];
            matrix_x(2,inode) = u[OverlapMap->LID(3*node+2)] + Mesh->nodes_coord[3*node+2];
            for (int iddl=0; iddl<3; ++iddl){
                Re(3*inode+iddl) = 0.0;
                Indices_tetra[3*inode+iddl] = 3*node+iddl;
                for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                    for (int jddl=0; jddl<3; ++jddl){
                        Ke(3*inode+iddl,3*jnode+jddl) = 0.0;
                    }
                }
            }
        }
                
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight = Mesh->gauss_weight_cells(gp);
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                dx_shape_functions(inode,0) = Mesh->DX_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,1) = Mesh->DY_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,2) = Mesh->DZ_N_tetra(gp+n_gauss_points*inode,e_lid);
            }
            
            error=deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
            compute_B_matrices(deformation_gradient,dx_shape_functions,matrix_B,matrix_BG);
            
            get_material_parameters(e_lid, gp);
            get_constitutive_tensors(deformation_gradient, scd_piola_stress, tangent_matrix);
            
            for (unsigned int i=0; i<3; ++i){
                block_scd_piola_stress(3*i+0,3*i+0) = scd_piola_stress(0);
                block_scd_piola_stress(3*i+0,3*i+1) = scd_piola_stress(5);
                block_scd_piola_stress(3*i+0,3*i+2) = scd_piola_stress(4);
                block_scd_piola_stress(3*i+1,3*i+0) = scd_piola_stress(5);
                block_scd_piola_stress(3*i+1,3*i+1) = scd_piola_stress(1);
                block_scd_piola_stress(3*i+1,3*i+2) = scd_piola_stress(3);
                block_scd_piola_stress(3*i+2,3*i+0) = scd_piola_stress(4);
                block_scd_piola_stress(3*i+2,3*i+1) = scd_piola_stress(3);
                block_scd_piola_stress(3*i+2,3*i+2) = scd_piola_stress(2);
            }
            
            error=Re.Multiply('T','N',-gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,scd_piola_stress,1.0);
            
            error=B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,tangent_matrix,0.0);
            error=Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
            
            error=BG_times_BSPS.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_BG,block_scd_piola_stress,0.0);
            error=Ke.Multiply('N','N',1.0,BG_times_BSPS,matrix_BG,1.0);
        }
        
        for (unsigned int i=0; i<3*Mesh->el_type; ++i){
            error=F.SumIntoGlobalValues(1, &Indices_tetra[i], &Re(i));
            for (unsigned int j=0; j<3*Mesh->el_type; ++j){
                error=K.SumIntoGlobalValues(1, &Indices_tetra[i], 1, &Indices_tetra[j], &Ke(i,j));
            }
        }
    }
    
}

void hyperelasticity_setup::material_stiffness_and_rhs_static_condensation(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    
    int node, e_gid, catch_error;
    int n_gauss_points = Mesh->n_gauss_cells;
    double gauss_weight;
    double det, theta, pressure, dpressure, coeff_kp;
    
    int *Indices_tetra;
    Indices_tetra = new int [3*Mesh->el_type];
    
    Epetra_SerialDenseMatrix Ke(3*Mesh->el_type,3*Mesh->el_type), Ke_vol(3*Mesh->el_type,3*Mesh->el_type);
    Epetra_SerialDenseMatrix Kg(3*Mesh->el_type,3*Mesh->el_type), Kg_vol(3*Mesh->el_type,3*Mesh->el_type);
    Epetra_SerialDenseMatrix Kp(3*Mesh->el_type,3*Mesh->el_type);
    Epetra_SerialDenseVector Re(3*Mesh->el_type), Re_vol(3*Mesh->el_type), Re_vol2(3*Mesh->el_type);
    
    Epetra_SerialDenseMatrix deformation_gradient(3,3);
    Epetra_SerialDenseMatrix tangent_matrix_isc(6,6);
    Epetra_SerialDenseMatrix tangent_matrix_vol(6,6);
    Epetra_SerialDenseVector scd_piola_stress_isc(6);
    Epetra_SerialDenseVector scd_piola_stress_vol(6);
    Epetra_SerialDenseVector piola_vol(6);
    Epetra_SerialDenseVector inverse_cauchy(6);
    Epetra_SerialDenseMatrix block_scd_piola_stress_isc(9,9);
    Epetra_SerialDenseMatrix block_scd_piola_stress_vol(9,9);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
    
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
    Epetra_SerialDenseMatrix B_times_TM(3*Mesh->el_type,6);
    Epetra_SerialDenseMatrix matrix_BG(9,3*Mesh->el_type);
    Epetra_SerialDenseMatrix BG_times_BSPS(3*Mesh->el_type,9);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            matrix_x(0,inode) = u[OverlapMap->LID(3*node+0)] + Mesh->nodes_coord[3*node+0];
            matrix_x(1,inode) = u[OverlapMap->LID(3*node+1)] + Mesh->nodes_coord[3*node+1];
            matrix_x(2,inode) = u[OverlapMap->LID(3*node+2)] + Mesh->nodes_coord[3*node+2];
            
            for (int iddl=0; iddl<3; ++iddl){
                Re(3*inode+iddl) = 0.0;
                Re_vol(3*inode+iddl) = 0.0;
                Re_vol2(3*inode+iddl) = 0.0;
                Indices_tetra[3*inode+iddl] = 3*node+iddl;
                for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                    for (int jddl=0; jddl<3; ++jddl){
                        Ke(3*inode+iddl,3*jnode+jddl) = 0.0;
                        Ke_vol(3*inode+iddl,3*jnode+jddl) = 0.0;
                        Kg(3*inode+iddl,3*jnode+jddl) = 0.0;
                        Kg_vol(3*inode+iddl,3*jnode+jddl) = 0.0;
                    }
                }
            }
        }
        
        theta = 0.0;
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight = Mesh->gauss_weight_cells(gp);
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                dx_shape_functions(inode,0) = Mesh->DX_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,1) = Mesh->DY_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,2) = Mesh->DZ_N_tetra(gp+n_gauss_points*inode,e_lid);
            }
            
            deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
            compute_B_matrices(deformation_gradient,dx_shape_functions,matrix_B,matrix_BG);
            
            get_material_parameters(e_lid, gp);
            get_constitutive_tensors_static_condensation(deformation_gradient, det, inverse_cauchy, scd_piola_stress_isc, scd_piola_stress_vol, tangent_matrix_isc, tangent_matrix_vol);
            
            for (unsigned int i=0; i<3; ++i){
                block_scd_piola_stress_isc(3*i+0,3*i+0) = scd_piola_stress_isc(0);
                block_scd_piola_stress_isc(3*i+0,3*i+1) = scd_piola_stress_isc(5);
                block_scd_piola_stress_isc(3*i+0,3*i+2) = scd_piola_stress_isc(4);
                block_scd_piola_stress_isc(3*i+1,3*i+0) = scd_piola_stress_isc(5);
                block_scd_piola_stress_isc(3*i+1,3*i+1) = scd_piola_stress_isc(1);
                block_scd_piola_stress_isc(3*i+1,3*i+2) = scd_piola_stress_isc(3);
                block_scd_piola_stress_isc(3*i+2,3*i+0) = scd_piola_stress_isc(4);
                block_scd_piola_stress_isc(3*i+2,3*i+1) = scd_piola_stress_isc(3);
                block_scd_piola_stress_isc(3*i+2,3*i+2) = scd_piola_stress_isc(2);

                block_scd_piola_stress_vol(3*i+0,3*i+0) = scd_piola_stress_vol(0);
                block_scd_piola_stress_vol(3*i+0,3*i+1) = scd_piola_stress_vol(5);
                block_scd_piola_stress_vol(3*i+0,3*i+2) = scd_piola_stress_vol(4);
                block_scd_piola_stress_vol(3*i+1,3*i+0) = scd_piola_stress_vol(5);
                block_scd_piola_stress_vol(3*i+1,3*i+1) = scd_piola_stress_vol(1);
                block_scd_piola_stress_vol(3*i+1,3*i+2) = scd_piola_stress_vol(3);
                block_scd_piola_stress_vol(3*i+2,3*i+0) = scd_piola_stress_vol(4);
                block_scd_piola_stress_vol(3*i+2,3*i+1) = scd_piola_stress_vol(3);
                block_scd_piola_stress_vol(3*i+2,3*i+2) = scd_piola_stress_vol(2);
            }
            
            for (unsigned int j=0; j<6; ++j){
                piola_vol(j) = det*inverse_cauchy(j);
            }
            
            theta += gauss_weight*det*Mesh->detJac_tetra(e_lid,gp);
                        
            Re.Multiply('T','N',-gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,scd_piola_stress_isc,1.0);
            Re_vol.Multiply('T','N',-gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,scd_piola_stress_vol,1.0);
            Re_vol2.Multiply('T','N',-gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,piola_vol,1.0);
            
            B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,tangent_matrix_isc,0.0);
            Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
            B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,tangent_matrix_vol,0.0);
            Ke_vol.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
            
            BG_times_BSPS.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_BG,block_scd_piola_stress_isc,0.0);
            Kg.Multiply('N','N',1.0,BG_times_BSPS,matrix_BG,1.0);
            BG_times_BSPS.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_BG,block_scd_piola_stress_vol,0.0);
            Kg_vol.Multiply('N','N',1.0,BG_times_BSPS,matrix_BG,1.0);
        }
        
        theta = theta/Mesh->vol_tetra(e_lid);
        get_internal_pressure(theta,pressure,dpressure);
        coeff_kp = dpressure/Mesh->vol_tetra(e_lid);
        
        Kp.Multiply('N','T',coeff_kp,Re_vol,Re_vol2,0.0);
        Re_vol.Scale(pressure);
        Ke_vol.Scale(pressure);
        Kg_vol.Scale(pressure);
        
        Re += Re_vol;
        Ke += Ke_vol;
        Ke += Kg;
        Ke += Kg_vol;
        Ke += Kp;
                
        for (unsigned int i=0; i<3*Mesh->el_type; ++i){
            F.SumIntoGlobalValues(1, &Indices_tetra[i], &Re(i));
            for (unsigned int j=0; j<3*Mesh->el_type; ++j){
                catch_error = K.SumIntoGlobalValues(1, &Indices_tetra[i], 1, &Indices_tetra[j], &Ke(i,j));
            }
        }
        
    }
    
    /*bool callFillComplete = false;
    K.GlobalAssemble(callFillComplete);
    F.GlobalAssemble();*/
}

void hyperelasticity_setup::force_stiffness_rhs_live_pressure(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    
    double kf;
    int n_gauss_points = Mesh->n_gauss_faces;
    double gauss_weight; // = Mesh->gauss_weight_tri;
    
    unsigned int e_gid;
    int node;
    int * Indices_tri;
    Indices_tri = new int [3*Mesh->face_type];
    
    Epetra_SerialDenseMatrix d_shape_functions(Mesh->face_type,2);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->face_type);
    Epetra_SerialDenseMatrix dxi_matrix_x(3,2);
    Epetra_SerialDenseMatrix k1(Mesh->face_type,Mesh->face_type), k2(Mesh->face_type,Mesh->face_type), k3(Mesh->face_type,Mesh->face_type);
    Epetra_SerialDenseVector force(3*Mesh->face_type);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_faces; ++e_lid){
        e_gid = Mesh->local_faces[e_lid];
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            node = Mesh->faces_nodes[Mesh->face_type*e_gid+inode];
            Indices_tri[3*inode+0] = 3*node+0;
            Indices_tri[3*inode+1] = 3*node+1;
            Indices_tri[3*inode+2] = 3*node+2;
            matrix_x(0,inode) = Mesh->nodes_coord[3*node+0] + u[OverlapMap->LID(3*node+0)];
            matrix_x(1,inode) = Mesh->nodes_coord[3*node+1] + u[OverlapMap->LID(3*node+1)];
            matrix_x(2,inode) = Mesh->nodes_coord[3*node+2] + u[OverlapMap->LID(3*node+2)];
            force(3*inode+0) = 0.0;
            force(3*inode+1) = 0.0;
            force(3*inode+2) = 0.0;
            for (unsigned int jnode=0; jnode<Mesh->face_type; ++jnode){
                k1(inode,jnode) = 0.0;
                k2(inode,jnode) = 0.0;
                k3(inode,jnode) = 0.0;
            }
        }
        
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight = Mesh->gauss_weight_faces(gp);
            for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
                d_shape_functions(inode,0) = Mesh->D1_N_tri(gp,inode);
                d_shape_functions(inode,1) = Mesh->D2_N_tri(gp,inode);
            }
            dxi_matrix_x.Multiply('N','N',1.0,matrix_x,d_shape_functions,0.0);
            
            for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
                force(3*inode+0) += gauss_weight*pressure_load*Mesh->N_tri(gp,inode)*(dxi_matrix_x(1,0)*dxi_matrix_x(2,1) - dxi_matrix_x(2,0)*dxi_matrix_x(1,1));
                force(3*inode+1) += gauss_weight*pressure_load*Mesh->N_tri(gp,inode)*(dxi_matrix_x(2,0)*dxi_matrix_x(0,1) - dxi_matrix_x(0,0)*dxi_matrix_x(2,1));
                force(3*inode+2) += gauss_weight*pressure_load*Mesh->N_tri(gp,inode)*(dxi_matrix_x(0,0)*dxi_matrix_x(1,1) - dxi_matrix_x(1,0)*dxi_matrix_x(0,1));
                
                for (unsigned int jnode=0; jnode<Mesh->face_type; ++jnode){
                    k1(inode,jnode) += 0.5*gauss_weight*pressure_load*( dxi_matrix_x(0,0)*(d_shape_functions(inode,1)*Mesh->N_tri(gp,jnode) - d_shape_functions(jnode,1)*Mesh->N_tri(gp,inode)) - dxi_matrix_x(0,1)*(d_shape_functions(inode,0)*Mesh->N_tri(gp,jnode) - d_shape_functions(jnode,0)*Mesh->N_tri(gp,inode)) );
                    k2(inode,jnode) += 0.5*gauss_weight*pressure_load*( dxi_matrix_x(1,0)*(d_shape_functions(inode,1)*Mesh->N_tri(gp,jnode) - d_shape_functions(jnode,1)*Mesh->N_tri(gp,inode)) - dxi_matrix_x(1,1)*(d_shape_functions(inode,0)*Mesh->N_tri(gp,jnode) - d_shape_functions(jnode,0)*Mesh->N_tri(gp,inode)) );
                    k3(inode,jnode) += 0.5*gauss_weight*pressure_load*( dxi_matrix_x(2,0)*(d_shape_functions(inode,1)*Mesh->N_tri(gp,jnode) - d_shape_functions(jnode,1)*Mesh->N_tri(gp,inode)) - dxi_matrix_x(2,1)*(d_shape_functions(inode,0)*Mesh->N_tri(gp,jnode) - d_shape_functions(jnode,0)*Mesh->N_tri(gp,inode)) );
                    //General case (leads to a non symmetric stiffness)
                    //k1(inode,jnode) += gauss_weight*pressure_load*Mesh->N_tri(gp,inode)*(d_shape_functions(jnode,0)*dxi_matrix_x(0,1)-d_shape_functions(jnode,1)*dxi_matrix_x(0,0));
                    //k2(inode,jnode) += gauss_Weight*pressure_load*Mesh->N_tri(gp,inode)*(d_shape_functions(jnode,0)*dxi_matrix_x(1,1)-d_shape_functions(jnode,1)*dxi_matrix_x(1,0));
                    //k3(inode,jnode) += gauss_weight*pressure_load*Mesh->N_tri(gp,inode)*(d_shape_functions(jnode,0)*dxi_matrix_x(2,1)-d_shape_functions(jnode,1)*dxi_matrix_x(2,0));
                }
            }
        }
        
        for (unsigned int i=0; i<3*Mesh->face_type; ++i){
            F.SumIntoGlobalValues(1, &Indices_tri[i], &force(i));
        }
        
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            for (unsigned int jnode=0; jnode<Mesh->face_type; ++jnode){
                kf = k1(inode,jnode);
                K.SumIntoGlobalValues(1, &Indices_tri[3*inode+2], 1, &Indices_tri[3*jnode+1], &kf);
                kf = -kf;
                K.SumIntoGlobalValues(1, &Indices_tri[3*inode+1], 1, &Indices_tri[3*jnode+2], &kf);
                
                kf = k2(inode,jnode);
                K.SumIntoGlobalValues(1, &Indices_tri[3*inode+0], 1, &Indices_tri[3*jnode+2], &kf);
                kf = -kf;
                K.SumIntoGlobalValues(1, &Indices_tri[3*inode+2], 1, &Indices_tri[3*jnode+0], &kf);
                
                kf = k3(inode,jnode);
                K.SumIntoGlobalValues(1, &Indices_tri[3*inode+1], 1, &Indices_tri[3*jnode+0], &kf);
                kf = -kf;
                K.SumIntoGlobalValues(1, &Indices_tri[3*inode+0], 1, &Indices_tri[3*jnode+1], &kf);
            }
        }
        
    }
    
    /*bool callFillComplete = false;
    K.GlobalAssemble(callFillComplete);
    F.GlobalAssemble();*/
}

void hyperelasticity_setup::force_dead_pressure(Epetra_FEVector & F){

    int node;
    int* Indices_tri;
    unsigned int e_gid;
    Indices_tri = new int [3*Mesh->face_type];
    
    int n_gauss_points = Mesh->n_gauss_faces;
    double gauss_weight;// = Mesh->gauss_weight_tri;
    
    Epetra_SerialDenseVector dead_pressure_load(3);
    Epetra_SerialDenseVector force(3*Mesh->face_type);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_faces; ++e_lid){
        e_gid  = Mesh->local_faces[e_lid];
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            node = Mesh->faces_nodes[Mesh->face_type*e_gid+inode];
            Indices_tri[3*inode] = 3*node;
            Indices_tri[3*inode+1] = 3*node+1;
            Indices_tri[3*inode+2] = 3*node+2;
            for (unsigned int iddl=0; iddl<3; ++iddl){
                force(3*inode+iddl) = 0.0;
            }
        }
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight = Mesh->gauss_weight_faces(gp);
            for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
                for (unsigned int iddl=0; iddl<3; ++iddl){
                    force(3*inode+iddl) += gauss_weight*dead_pressure_load(iddl)*Mesh->N_tri(gp,inode)*Mesh->detJac_tri(e_lid,gp);
                }
            }
        }
        
        for (unsigned int inode=0; inode<6; ++inode){
            for (unsigned int iddl=0; iddl<3; ++iddl){
                F.SumIntoGlobalValues(1, &Indices_tri[3*inode+iddl], &force(3*inode+iddl));
            }
        }
        
    }
}

void hyperelasticity_setup::compute_green_lagrange(Epetra_Vector & x, double & xi, double & eta, double & zeta){
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    Epetra_Map CellsMap(-1,Mesh->n_local_cells,&Mesh->local_cells[0],0,*Comm);
    Epetra_Vector green_lagrange(CellsMap);
    
    int node, e_gid;
    double det_jac_tetra;
        
    Epetra_SerialDenseMatrix deformation_gradient(3,3);
    Epetra_SerialDenseMatrix right_cauchy(3,3);
    Epetra_SerialDenseMatrix matrix_X(3,Mesh->el_type);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix D(Mesh->el_type,3), DX(Mesh->el_type,3);
    Epetra_SerialDenseMatrix JacobianMatrix(3,3);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            matrix_X(0,inode) = Mesh->nodes_coord[3*node+0];
            matrix_X(1,inode) = Mesh->nodes_coord[3*node+1];
            matrix_X(2,inode) = Mesh->nodes_coord[3*node+2];
            matrix_x(0,inode) = u[OverlapMap->LID(3*node+0)] + Mesh->nodes_coord[3*node+0];
            matrix_x(1,inode) = u[OverlapMap->LID(3*node+1)] + Mesh->nodes_coord[3*node+1];
            matrix_x(2,inode) = u[OverlapMap->LID(3*node+2)] + Mesh->nodes_coord[3*node+2];
        }
        
        switch (Mesh->el_type){
            case 4:
                tetra4::d_shape_functions(D, xi, eta, zeta);
                break;
            case 8:
                hexa8::d_shape_functions(D, xi, eta, zeta);
                break;
            case 10:
                tetra10::d_shape_functions(D, xi, eta, zeta);
                break;
        }
        jacobian_matrix(matrix_X,D,JacobianMatrix);
        jacobian_det(JacobianMatrix,det_jac_tetra);
        dX_shape_functions(D,JacobianMatrix,det_jac_tetra,dx_shape_functions);
        
        deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
        right_cauchy.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);
        green_lagrange[e_lid] = 0.5*(right_cauchy(1,1)-1.0);
    }
    
    int NumTargetElements = 0;
    if (Comm->MyPID()==0){
        NumTargetElements = Mesh->n_cells;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export ExportOnRoot(CellsMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(green_lagrange,ExportOnRoot,Insert);
    
    int error = EpetraExt::MultiVectorToMatrixMarketFile("green_lagrange.mtx",lhs_root,0,0,false);
}

void hyperelasticity_setup::compute_cauchy_stress(Epetra_Vector & x, std::string & filename){
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    int node, e_gid;
    int n_gauss_points = Mesh->n_gauss_cells;
    double det_jac_tetra, det, gauss_weight;
    
    Epetra_SerialDenseMatrix deformation_gradient(3,3);
    Epetra_SerialDenseMatrix right_cauchy(3,3);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix piola_stress(3,3);
    Epetra_SerialDenseMatrix cauchy_stress(3,3);
    Epetra_SerialDenseMatrix dg_times_ps(3,3);
    
    std::vector<int> local_gauss_points;
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        for (unsigned int j=0; j<Mesh->n_gauss_cells; ++j){
            local_gauss_points.push_back(Mesh->n_gauss_cells*e_gid+j);
        }
    }
    Epetra_Map CellsMap(-1,Mesh->n_local_cells,&local_gauss_points[0],0,*Comm);
    Epetra_Vector sigma11(CellsMap);
    Epetra_Vector sigma22(CellsMap);
    Epetra_Vector sigma33(CellsMap);
    Epetra_Vector sigma12(CellsMap);
    Epetra_Vector sigma13(CellsMap);
    Epetra_Vector sigma23(CellsMap);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            matrix_x(0,inode) = u[OverlapMap->LID(3*node+0)] + Mesh->nodes_coord[3*node+0];
            matrix_x(1,inode) = u[OverlapMap->LID(3*node+1)] + Mesh->nodes_coord[3*node+1];
            matrix_x(2,inode) = u[OverlapMap->LID(3*node+2)] + Mesh->nodes_coord[3*node+2];
        }
        
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight = Mesh->gauss_weight_cells(gp);
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                dx_shape_functions(inode,0) = Mesh->DX_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,1) = Mesh->DY_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,2) = Mesh->DZ_N_tetra(gp+n_gauss_points*inode,e_lid);
            }
            
            deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
            
            get_material_parameters(e_lid, gp);
            get_stress_for_recover(deformation_gradient, det, piola_stress);
            dg_times_ps.Multiply('N','N',1.0,deformation_gradient,piola_stress,0.0);
            cauchy_stress.Multiply('N','T',1.0,dg_times_ps,deformation_gradient,0.0);
            
            sigma11[Mesh->n_gauss_cells*e_lid+gp] = cauchy_stress(0,0);
            sigma22[Mesh->n_gauss_cells*e_lid+gp] = cauchy_stress(1,1);
            sigma33[Mesh->n_gauss_cells*e_lid+gp] = cauchy_stress(2,2);
            sigma12[Mesh->n_gauss_cells*e_lid+gp] = cauchy_stress(0,1);
            sigma13[Mesh->n_gauss_cells*e_lid+gp] = cauchy_stress(0,2);
            sigma23[Mesh->n_gauss_cells*e_lid+gp] = cauchy_stress(1,2);
        }
    }
    
    int NumTargetElements = 0;
    if (Comm->MyPID()==0){
        NumTargetElements = Mesh->n_cells;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export ExportOnRoot(CellsMap,MapOnRoot);
    
    Epetra_MultiVector lhs_root11(MapOnRoot,true);
    lhs_root11.Export(sigma11,ExportOnRoot,Insert);
    std::string file11 = "sigma11" + filename;
    int error11 = EpetraExt::MultiVectorToMatrixMarketFile(file11.c_str(),lhs_root11,0,0,false);
    
    Epetra_MultiVector lhs_root22(MapOnRoot,true);
    lhs_root22.Export(sigma22,ExportOnRoot,Insert);
    std::string file22 = "sigma22" + filename;
    int error22 = EpetraExt::MultiVectorToMatrixMarketFile(file22.c_str(),lhs_root22,0,0,false);
    
    Epetra_MultiVector lhs_root33(MapOnRoot,true);
    lhs_root33.Export(sigma33,ExportOnRoot,Insert);
    std::string file33 = "sigma33" + filename;
    int error33 = EpetraExt::MultiVectorToMatrixMarketFile(file33.c_str(),lhs_root33,0,0,false);
    
    Epetra_MultiVector lhs_root12(MapOnRoot,true);
    lhs_root12.Export(sigma12,ExportOnRoot,Insert);
    std::string file12 = "sigma12" + filename;
    int error12 = EpetraExt::MultiVectorToMatrixMarketFile(file12.c_str(),lhs_root12,0,0,false);
    
    Epetra_MultiVector lhs_root13(MapOnRoot,true);
    lhs_root13.Export(sigma13,ExportOnRoot,Insert);
    std::string file13 = "sigma13" + filename;
    int error13 = EpetraExt::MultiVectorToMatrixMarketFile(file13.c_str(),lhs_root13,0,0,false);
    
    Epetra_MultiVector lhs_root23(MapOnRoot,true);
    lhs_root23.Export(sigma23,ExportOnRoot,Insert);
    std::string file23 = "sigma23" + filename;
    int error23 = EpetraExt::MultiVectorToMatrixMarketFile(file23.c_str(),lhs_root23,0,0,false);
    
}

void hyperelasticity_setup::compute_mean_cauchy_stress(Epetra_Vector & x, std::string & filename){
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    Epetra_Map CellsMap(-1,Mesh->n_local_cells,&Mesh->local_cells[0],0,*Comm);
    Epetra_Vector sigma11(CellsMap); sigma11.PutScalar(0.0);
    Epetra_Vector sigma22(CellsMap); sigma22.PutScalar(0.0);
    Epetra_Vector sigma33(CellsMap); sigma33.PutScalar(0.0);
    Epetra_Vector sigma12(CellsMap); sigma12.PutScalar(0.0);
    Epetra_Vector sigma13(CellsMap); sigma13.PutScalar(0.0);
    Epetra_Vector sigma23(CellsMap); sigma23.PutScalar(0.0);
    
    int node, e_gid;
    int n_gauss_points = Mesh->n_gauss_cells;
    double det_jac_tetra, det, gauss_weight, theta;
    
    Epetra_SerialDenseMatrix deformation_gradient(3,3);
    Epetra_SerialDenseMatrix right_cauchy(3,3);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix piola_stress(3,3);
    Epetra_SerialDenseMatrix cauchy_stress(3,3);
    Epetra_SerialDenseMatrix dg_times_ps(3,3);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            matrix_x(0,inode) = u[OverlapMap->LID(3*node+0)] + Mesh->nodes_coord[3*node+0];
            matrix_x(1,inode) = u[OverlapMap->LID(3*node+1)] + Mesh->nodes_coord[3*node+1];
            matrix_x(2,inode) = u[OverlapMap->LID(3*node+2)] + Mesh->nodes_coord[3*node+2];
        }
        
        theta = 0.0;
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight = Mesh->gauss_weight_cells(gp);
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                dx_shape_functions(inode,0) = Mesh->DX_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,1) = Mesh->DY_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,2) = Mesh->DZ_N_tetra(gp+n_gauss_points*inode,e_lid);
            }
        
            deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
            get_material_parameters(e_lid, gp);
            get_stress_for_recover(deformation_gradient, det, piola_stress);
            dg_times_ps.Multiply('N','N',1.0,deformation_gradient,piola_stress,0.0);
            cauchy_stress.Multiply('N','T',1.0,dg_times_ps,deformation_gradient,0.0);
            
            sigma11[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(0,0);
            sigma22[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(1,1);
            sigma33[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(2,2);
            sigma12[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(0,1);
            sigma13[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(0,2);
            sigma23[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(1,2);
            
            theta += gauss_weight*det*Mesh->detJac_tetra(e_lid,gp);
            
        }
        sigma11[e_lid] = sigma11[e_lid]/theta;
        sigma22[e_lid] = sigma22[e_lid]/theta;
        sigma33[e_lid] = sigma33[e_lid]/theta;
        sigma12[e_lid] = sigma12[e_lid]/theta;
        sigma13[e_lid] = sigma13[e_lid]/theta;
        sigma23[e_lid] = sigma23[e_lid]/theta;
    }
    
    int NumTargetElements = 0;
    if (Comm->MyPID()==0){
        NumTargetElements = Mesh->n_cells;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export ExportOnRoot(CellsMap,MapOnRoot);
    
    Epetra_MultiVector lhs_root11(MapOnRoot,true);
    lhs_root11.Export(sigma11,ExportOnRoot,Insert);
    std::string file11 = filename + "_sigma11.mtx";
    int error11 = EpetraExt::MultiVectorToMatrixMarketFile(file11.c_str(),lhs_root11,0,0,false);

    Epetra_MultiVector lhs_root22(MapOnRoot,true);
    lhs_root22.Export(sigma22,ExportOnRoot,Insert);
    std::string file22 = filename + "_sigma22.mtx";
    int error22 = EpetraExt::MultiVectorToMatrixMarketFile(file22.c_str(),lhs_root22,0,0,false);
    
    Epetra_MultiVector lhs_root33(MapOnRoot,true);
    lhs_root33.Export(sigma33,ExportOnRoot,Insert);
    std::string file33 = filename + "_sigma33.mtx";
    int error33 = EpetraExt::MultiVectorToMatrixMarketFile(file33.c_str(),lhs_root33,0,0,false);
    
    Epetra_MultiVector lhs_root12(MapOnRoot,true);
    lhs_root12.Export(sigma12,ExportOnRoot,Insert);
    std::string file12 = filename + "_sigma12.mtx";
    int error12 = EpetraExt::MultiVectorToMatrixMarketFile(file12.c_str(),lhs_root12,0,0,false);
    
    Epetra_MultiVector lhs_root13(MapOnRoot,true);
    lhs_root13.Export(sigma13,ExportOnRoot,Insert);
    std::string file13 = filename + "_sigma13.mtx";
    int error13 = EpetraExt::MultiVectorToMatrixMarketFile(file13.c_str(),lhs_root13,0,0,false);
    
    Epetra_MultiVector lhs_root23(MapOnRoot,true);
    lhs_root23.Export(sigma23,ExportOnRoot,Insert);
    std::string file23 = filename + "_sigma23.mtx";
    int error23 = EpetraExt::MultiVectorToMatrixMarketFile(file23.c_str(),lhs_root23,0,0,false);
}

void hyperelasticity_setup::compute_mean_von_mises_stress(Epetra_Vector & x, std::string & filename){
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    Epetra_Map CellsMap(-1,Mesh->n_local_cells,&Mesh->local_cells[0],0,*Comm);
    Epetra_Vector von_mises(CellsMap);
    
    int node, e_gid;
    int n_gauss_points = Mesh->n_gauss_cells;
    double det_jac_tetra, det, gauss_weight, theta;
    
    Epetra_SerialDenseMatrix deformation_gradient(3,3);
    Epetra_SerialDenseMatrix right_cauchy(3,3);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix piola_stress(3,3);
    Epetra_SerialDenseMatrix cauchy_stress(3,3);
    Epetra_SerialDenseMatrix dg_times_ps(3,3);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            matrix_x(0,inode) = u[OverlapMap->LID(3*node+0)] + Mesh->nodes_coord[3*node+0];
            matrix_x(1,inode) = u[OverlapMap->LID(3*node+1)] + Mesh->nodes_coord[3*node+1];
            matrix_x(2,inode) = u[OverlapMap->LID(3*node+2)] + Mesh->nodes_coord[3*node+2];
        }
        
        theta = 0.0;
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight = Mesh->gauss_weight_cells(gp);
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                dx_shape_functions(inode,0) = Mesh->DX_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,1) = Mesh->DY_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,2) = Mesh->DZ_N_tetra(gp+n_gauss_points*inode,e_lid);
            }
            
            deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
            
            get_material_parameters(e_lid, gp);
            get_stress_for_recover(deformation_gradient, det, piola_stress);
            dg_times_ps.Multiply('N','N',1.0,deformation_gradient,piola_stress,0.0);
            cauchy_stress.Multiply('N','T',1.0,dg_times_ps,deformation_gradient,0.0);
            
            von_mises[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)* ( std::sqrt( (cauchy_stress(0,0)-cauchy_stress(1,1))*(cauchy_stress(0,0)-cauchy_stress(1,1)) + (cauchy_stress(1,1)-cauchy_stress(2,2))*(cauchy_stress(1,1)-cauchy_stress(2,2)) + (cauchy_stress(2,2)-cauchy_stress(0,0))*(cauchy_stress(2,2)-cauchy_stress(0,0)) + 6.0*(cauchy_stress(1,2)*cauchy_stress(1,2) + cauchy_stress(2,0)*cauchy_stress(2,0) + cauchy_stress(0,1)*cauchy_stress(0,1)) ) );
            theta += gauss_weight*det*Mesh->detJac_tetra(e_lid,gp);
            
        }
        von_mises[e_lid] = von_mises[e_lid]/theta;
    }
    
    int NumTargetElements = 0;
    if (Comm->MyPID()==0){
        NumTargetElements = Mesh->n_cells;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export ExportOnRoot(CellsMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(von_mises,ExportOnRoot,Insert);
    
    int error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    
}

void hyperelasticity_setup::compute_B_matrices(Epetra_SerialDenseMatrix & F, Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B, Epetra_SerialDenseMatrix & BG){
    
    for (unsigned inode=0; inode<Mesh->el_type; ++inode){
        B(0,3*inode) = F(0,0)*dx_shape_functions(inode,0);
        B(0,3*inode+1) = F(1,0)*dx_shape_functions(inode,0);
        B(0,3*inode+2) = F(2,0)*dx_shape_functions(inode,0);
    
        B(1,3*inode) = F(0,1)*dx_shape_functions(inode,1);
        B(1,3*inode+1) = F(1,1)*dx_shape_functions(inode,1);
        B(1,3*inode+2) = F(2,1)*dx_shape_functions(inode,1);
    
        B(2,3*inode) = F(0,2)*dx_shape_functions(inode,2);
        B(2,3*inode+1) = F(1,2)*dx_shape_functions(inode,2);
        B(2,3*inode+2) = F(2,2)*dx_shape_functions(inode,2);
    
        B(3,3*inode) = F(0,1)*dx_shape_functions(inode,2) + F(0,2)*dx_shape_functions(inode,1);
        B(3,3*inode+1) = F(1,1)*dx_shape_functions(inode,2) + F(1,2)*dx_shape_functions(inode,1);
        B(3,3*inode+2) = F(2,1)*dx_shape_functions(inode,2) + F(2,2)*dx_shape_functions(inode,1);
    
        B(4,3*inode) = F(0,0)*dx_shape_functions(inode,2) + F(0,2)*dx_shape_functions(inode,0);
        B(4,3*inode+1) = F(1,0)*dx_shape_functions(inode,2) + F(1,2)*dx_shape_functions(inode,0);
        B(4,3*inode+2) = F(2,0)*dx_shape_functions(inode,2) + F(2,2)*dx_shape_functions(inode,0);
    
        B(5,3*inode) = F(0,0)*dx_shape_functions(inode,1) + F(0,1)*dx_shape_functions(inode,0);
        B(5,3*inode+1) = F(1,0)*dx_shape_functions(inode,1) + F(1,1)*dx_shape_functions(inode,0);
        B(5,3*inode+2) = F(2,0)*dx_shape_functions(inode,1) + F(2,1)*dx_shape_functions(inode,0);
    
        BG(0,3*inode) = dx_shape_functions(inode,0);
        BG(1,3*inode) = dx_shape_functions(inode,1);
        BG(2,3*inode) = dx_shape_functions(inode,2);
    
        BG(3,3*inode+1) = dx_shape_functions(inode,0);
        BG(4,3*inode+1) = dx_shape_functions(inode,1);
        BG(5,3*inode+1) = dx_shape_functions(inode,2);
    
        BG(6,3*inode+2) = dx_shape_functions(inode,0);
        BG(7,3*inode+2) = dx_shape_functions(inode,1);
        BG(8,3*inode+2) = dx_shape_functions(inode,2);
    }
    
}

/*void hyperelasticity_setup::material_stiffness_and_rhs_static_condensation_twosteps(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
 
 int node, e_gid, catch_error;
 int n_gauss_points = Mesh->n_gauss_cells;
 double gauss_weight;
 double det, theta, pressure, dpressure, coeff_kp;
 
 int *Indices_tetra;
 Indices_tetra = new int [3*Mesh->el_type];
 
 Epetra_SerialDenseMatrix Ke(3*Mesh->el_type,3*Mesh->el_type);
 Epetra_SerialDenseMatrix Kg(3*Mesh->el_type,3*Mesh->el_type);
 Epetra_SerialDenseMatrix Kp(3*Mesh->el_type,3*Mesh->el_type);
 Epetra_SerialDenseVector Re(3*Mesh->el_type);
 
 Epetra_SerialDenseMatrix deformation_gradient(3,3);
 Epetra_SerialDenseMatrix tangent_matrix(6,6);
 Epetra_SerialDenseVector scd_piola_stress(6);
 Epetra_SerialDenseMatrix block_scd_piola_stress(9,9);
 Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
 
 Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
 Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
 Epetra_SerialDenseMatrix B_times_TM(3*Mesh->el_type,6);
 Epetra_SerialDenseMatrix matrix_BG(9,3*Mesh->el_type);
 Epetra_SerialDenseMatrix BG_times_BSPS(3*Mesh->el_type,9);
 
 for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
 e_gid = Mesh->local_cells[e_lid];
 
 for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
 node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
 matrix_x(0,inode) = u[OverlapMap->LID(3*node+0)] + Mesh->nodes_coord[3*node+0];
 matrix_x(1,inode) = u[OverlapMap->LID(3*node+1)] + Mesh->nodes_coord[3*node+1];
 matrix_x(2,inode) = u[OverlapMap->LID(3*node+2)] + Mesh->nodes_coord[3*node+2];
 
 for (int iddl=0; iddl<3; ++iddl){
 Re(3*inode+iddl) = 0.0;
 Indices_tetra[3*inode+iddl] = 3*node+iddl;
 for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
 for (int jddl=0; jddl<3; ++jddl){
 Ke(3*inode+iddl,3*jnode+jddl) = 0.0;
 Kg(3*inode+iddl,3*jnode+jddl) = 0.0;
 }
 }
 }
 }
 
 theta = 0.0;
 theta = theta/Mesh->vol_tetra(e_lid);
 coeff_kp = dpressure/Mesh->vol_tetra(e_lid);
 get_internal_pressure(theta,pressure,dpressure);
 
 for (unsigned int gp=0; gp<n_gauss_points; ++gp){
 gauss_weight = Mesh->gauss_weight_cells(gp);
 for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
 dx_shape_functions(inode,0) = Mesh->DX_N_tetra(gp+n_gauss_points*inode,e_lid);
 dx_shape_functions(inode,1) = Mesh->DY_N_tetra(gp+n_gauss_points*inode,e_lid);
 dx_shape_functions(inode,2) = Mesh->DZ_N_tetra(gp+n_gauss_points*inode,e_lid);
 }
 
 deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
 compute_B_matrices(deformation_gradient,dx_shape_functions,matrix_B,matrix_BG);
 
 get_material_parameters(e_lid, gp);
 get_constitutive_tensors_static_condensation(deformation_gradient, scd_piola_stress, tangent_matrix);
 
 for (unsigned int i=0; i<3; ++i){
 block_scd_piola_stress(3*i+0,3*i+0) = scd_piola_stress(0);
 block_scd_piola_stress(3*i+0,3*i+1) = scd_piola_stress(5);
 block_scd_piola_stress(3*i+0,3*i+2) = scd_piola_stress(4);
 block_scd_piola_stress(3*i+1,3*i+0) = scd_piola_stress(5);
 block_scd_piola_stress(3*i+1,3*i+1) = scd_piola_stress(1);
 block_scd_piola_stress(3*i+1,3*i+2) = scd_piola_stress(3);
 block_scd_piola_stress(3*i+2,3*i+0) = scd_piola_stress(4);
 block_scd_piola_stress(3*i+2,3*i+1) = scd_piola_stress(3);
 block_scd_piola_stress(3*i+2,3*i+2) = scd_piola_stress(2);
 }
 
 //theta += gauss_weight*det*Mesh->detJac_tetra(e_lid,gp);
 
 Re.Multiply('T','N',-gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,scd_piola_stress,1.0);
 
 B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,tangent_matrix_isc,0.0);
 Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
 
 BG_times_BSPS.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_BG,block_scd_piola_stress,0.0);
 Kg.Multiply('N','N',1.0,BG_times_BSPS,matrix_BG,1.0);
 }
 
 Kp.Multiply('N','T',coeff_kp,Re_vol,Re_vol,0.0);
 Ke += Kg;
 Ke += Kp;
 
 for (unsigned int i=0; i<3*Mesh->el_type; ++i){
 F.SumIntoGlobalValues(1, &Indices_tetra[i], &Re(i));
 for (unsigned int j=0; j<3*Mesh->el_type; ++j){
 catch_error = K.SumIntoGlobalValues(1, &Indices_tetra[i], 1, &Indices_tetra[j], &Ke(i,j));
 }
 }
 
 }
 
 }*/
