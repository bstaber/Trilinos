/*
Brian Staber (brian.staber@gmail.com)
*/

#include "nearlyIncompressibleHyperelasticity.hpp"

nearlyIncompressibleHyperelasticity::nearlyIncompressibleHyperelasticity(){
}

nearlyIncompressibleHyperelasticity::~nearlyIncompressibleHyperelasticity(){
}

void nearlyIncompressibleHyperelasticity::assemblePureDirichlet_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){

    F.PutScalar(0.0);
    K.PutScalar(0.0);

    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);

    stiffnessRhsMaterialContribution(u,K,F);

    Comm->Barrier();

    K.GlobalAssemble();
    K.FillComplete();
    F.GlobalAssemble();
}
void nearlyIncompressibleHyperelasticity::assembleMixedDirichletDeformationDependentNeumann_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){

    F.PutScalar(0.0);
    K.PutScalar(0.0);

    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);

    stiffnessRhsMaterialContribution(u,K,F);
    stiffnessRhsPressureContribution(u,K,F);

    Comm->Barrier();
    K.GlobalAssemble();
    K.FillComplete();
    F.GlobalAssemble();
}

void nearlyIncompressibleHyperelasticity::stiffnessRhsMaterialContribution(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    int node, e_gid, catch_error;
    int n_gauss_points = Mesh->n_gauss_cells;
    double gauss_weight;
    double det, theta, pressure, dpressure, coeff_kp;

    int *Indexes;
    Indexes = new int [3*Mesh->el_type];

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
                Indexes[3*inode+iddl] = 3*node+iddl;
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
                dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,e_lid);
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

            theta += gauss_weight*det*Mesh->detJac_cells(e_lid,gp);

            Re.Multiply('T','N',-gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,scd_piola_stress_isc,1.0);
            Re_vol.Multiply('T','N',-gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,scd_piola_stress_vol,1.0);
            Re_vol2.Multiply('T','N',-gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,piola_vol,1.0);

            B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,tangent_matrix_isc,0.0);
            Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
            B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,tangent_matrix_vol,0.0);
            Ke_vol.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);

            BG_times_BSPS.Multiply('T','N',gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_BG,block_scd_piola_stress_isc,0.0);
            Kg.Multiply('N','N',1.0,BG_times_BSPS,matrix_BG,1.0);
            BG_times_BSPS.Multiply('T','N',gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_BG,block_scd_piola_stress_vol,0.0);
            Kg_vol.Multiply('N','N',1.0,BG_times_BSPS,matrix_BG,1.0);
        }

        theta = theta/Mesh->vol_cells(e_lid);
        get_internal_pressure(theta,pressure,dpressure);
        coeff_kp = dpressure/Mesh->vol_cells(e_lid);

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
            F.SumIntoGlobalValues(1, &Indexes[i], &Re(i));
            for (unsigned int j=0; j<3*Mesh->el_type; ++j){
                catch_error = K.SumIntoGlobalValues(1, &Indexes[i], 1, &Indexes[j], &Ke(i,j));
            }
        }

    }
    delete[] Indexes;
}
void nearlyIncompressibleHyperelasticity::stiffnessRhsPressureContribution(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F){

    double kf;
    int n_gauss_points = Mesh->n_gauss_faces;
    double gauss_weight;

    unsigned int e_gid;
    int node;
    int * Indexes;
    Indexes = new int [3*Mesh->face_type];

    Epetra_SerialDenseMatrix d_shape_functions(Mesh->face_type,2);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->face_type);
    Epetra_SerialDenseMatrix dxi_matrix_x(3,2);
    Epetra_SerialDenseMatrix k1(Mesh->face_type,Mesh->face_type), k2(Mesh->face_type,Mesh->face_type), k3(Mesh->face_type,Mesh->face_type);
    Epetra_SerialDenseVector force(3*Mesh->face_type);

    for (unsigned int e_lid=0; e_lid<Mesh->n_local_faces; ++e_lid){
        e_gid = Mesh->local_faces[e_lid];
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            node = Mesh->faces_nodes[Mesh->face_type*e_gid+inode];
            Indexes[3*inode+0] = 3*node+0;
            Indexes[3*inode+1] = 3*node+1;
            Indexes[3*inode+2] = 3*node+2;
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
                d_shape_functions(inode,0) = Mesh->D1_N_faces(gp,inode);
                d_shape_functions(inode,1) = Mesh->D2_N_faces(gp,inode);
            }
            dxi_matrix_x.Multiply('N','N',1.0,matrix_x,d_shape_functions,0.0);

            for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
                force(3*inode+0) += gauss_weight*pressure_load*Mesh->N_faces(gp,inode)*(dxi_matrix_x(1,0)*dxi_matrix_x(2,1) - dxi_matrix_x(2,0)*dxi_matrix_x(1,1));
                force(3*inode+1) += gauss_weight*pressure_load*Mesh->N_faces(gp,inode)*(dxi_matrix_x(2,0)*dxi_matrix_x(0,1) - dxi_matrix_x(0,0)*dxi_matrix_x(2,1));
                force(3*inode+2) += gauss_weight*pressure_load*Mesh->N_faces(gp,inode)*(dxi_matrix_x(0,0)*dxi_matrix_x(1,1) - dxi_matrix_x(1,0)*dxi_matrix_x(0,1));

                for (unsigned int jnode=0; jnode<Mesh->face_type; ++jnode){
                    k1(inode,jnode) += 0.5*gauss_weight*pressure_load*( dxi_matrix_x(0,0)*(d_shape_functions(inode,1)*Mesh->N_faces(gp,jnode) - d_shape_functions(jnode,1)*Mesh->N_faces(gp,inode)) - dxi_matrix_x(0,1)*(d_shape_functions(inode,0)*Mesh->N_faces(gp,jnode) - d_shape_functions(jnode,0)*Mesh->N_faces(gp,inode)) );
                    k2(inode,jnode) += 0.5*gauss_weight*pressure_load*( dxi_matrix_x(1,0)*(d_shape_functions(inode,1)*Mesh->N_faces(gp,jnode) - d_shape_functions(jnode,1)*Mesh->N_faces(gp,inode)) - dxi_matrix_x(1,1)*(d_shape_functions(inode,0)*Mesh->N_faces(gp,jnode) - d_shape_functions(jnode,0)*Mesh->N_faces(gp,inode)) );
                    k3(inode,jnode) += 0.5*gauss_weight*pressure_load*( dxi_matrix_x(2,0)*(d_shape_functions(inode,1)*Mesh->N_faces(gp,jnode) - d_shape_functions(jnode,1)*Mesh->N_faces(gp,inode)) - dxi_matrix_x(2,1)*(d_shape_functions(inode,0)*Mesh->N_faces(gp,jnode) - d_shape_functions(jnode,0)*Mesh->N_faces(gp,inode)) );
                    //General case (leads to a non symmetric stiffness)
                    //k1(inode,jnode) += gauss_weight*pressure_load*Mesh->N_faces(gp,inode)*(d_shape_functions(jnode,0)*dxi_matrix_x(0,1)-d_shape_functions(jnode,1)*dxi_matrix_x(0,0));
                    //k2(inode,jnode) += gauss_Weight*pressure_load*Mesh->N_faces(gp,inode)*(d_shape_functions(jnode,0)*dxi_matrix_x(1,1)-d_shape_functions(jnode,1)*dxi_matrix_x(1,0));
                    //k3(inode,jnode) += gauss_weight*pressure_load*Mesh->N_faces(gp,inode)*(d_shape_functions(jnode,0)*dxi_matrix_x(2,1)-d_shape_functions(jnode,1)*dxi_matrix_x(2,0));
                }
            }
        }

        for (unsigned int i=0; i<3*Mesh->face_type; ++i){
            F.SumIntoGlobalValues(1, &Indexes[i], &force(i));
        }

        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            for (unsigned int jnode=0; jnode<Mesh->face_type; ++jnode){
                kf = k1(inode,jnode);
                K.SumIntoGlobalValues(1, &Indexes[3*inode+2], 1, &Indexes[3*jnode+1], &kf);
                kf = -kf;
                K.SumIntoGlobalValues(1, &Indexes[3*inode+1], 1, &Indexes[3*jnode+2], &kf);

                kf = k2(inode,jnode);
                K.SumIntoGlobalValues(1, &Indexes[3*inode+0], 1, &Indexes[3*jnode+2], &kf);
                kf = -kf;
                K.SumIntoGlobalValues(1, &Indexes[3*inode+2], 1, &Indexes[3*jnode+0], &kf);

                kf = k3(inode,jnode);
                K.SumIntoGlobalValues(1, &Indexes[3*inode+1], 1, &Indexes[3*jnode+0], &kf);
                kf = -kf;
                K.SumIntoGlobalValues(1, &Indexes[3*inode+0], 1, &Indexes[3*jnode+1], &kf);
            }
        }

    }
    delete[] Indexes;
}
