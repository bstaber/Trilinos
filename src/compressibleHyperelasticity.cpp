#include "compressibleHyperelasticity.hpp"

compressibleHyperelasticity::compressibleHyperelasticity(){
}

compressibleHyperelasticity::~compressibleHyperelasticity(){
}

void compressibleHyperelasticity::assemblePureDirichlet_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    int error;
    F.PutScalar(0.0);
    K.PutScalar(0.0);
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    stiffnessRhs_homogeneousForcing(u,K,F);
    
    Comm->Barrier();
    error = K.GlobalAssemble();
    error = K.FillComplete();
    error = F.GlobalAssemble();
}
void compressibleHyperelasticity::assemblePureDirichlet_inhomogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    int error;
    F.PutScalar(0.0);
    K.PutScalar(0.0);
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    stiffnessRhs_inhomogeneousForcing(u,K,F);
    
    Comm->Barrier();
    error = K.GlobalAssemble();
    error = K.FillComplete();
    error = F.GlobalAssemble();
}
void compressibleHyperelasticity::assembleMixedDirichletNeumann_homogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    int error;
    F.PutScalar(0.0);
    K.PutScalar(0.0);
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    stiffnessRhs_homogeneousForcing(u,K,F);
    rhs_NeumannBoundaryCondition(F);
    
    Comm->Barrier();
    error = K.GlobalAssemble();
    error = K.FillComplete();
    error = F.GlobalAssemble();
}
void compressibleHyperelasticity::assembleMixedDirichletNeumann_inhomogeneousForcing(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    int error;
    F.PutScalar(0.0);
    K.PutScalar(0.0);
    
    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);
    
    stiffnessRhs_inhomogeneousForcing(u,K,F);
    rhs_NeumannBoundaryCondition(F);
    
    Comm->Barrier();
    error = K.GlobalAssemble();
    error = K.FillComplete();
    error = F.GlobalAssemble();
}

void compressibleHyperelasticity::stiffnessRhs_homogeneousForcing(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    int node, e_gid, error;
    int n_gauss_points = Mesh->n_gauss_cells;
    double gauss_weight;
    
    int *Indexes;
    Indexes = new int [3*Mesh->el_type];
    
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
                Indexes[3*inode+iddl] = 3*node+iddl;
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
            
            error = Re.Multiply('T','N',-gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,scd_piola_stress,1.0);
            
            error = B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,tangent_matrix,0.0);
            error = Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
            
            error = BG_times_BSPS.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_BG,block_scd_piola_stress,0.0);
            error = Ke.Multiply('N','N',1.0,BG_times_BSPS,matrix_BG,1.0);
        }
        
        for (unsigned int i=0; i<3*Mesh->el_type; ++i){
            error = F.SumIntoGlobalValues(1, &Indexes[i], &Re(i));
            for (unsigned int j=0; j<3*Mesh->el_type; ++j){
                error = K.SumIntoGlobalValues(1, &Indexes[i], 1, &Indexes[j], &Ke(i,j));
            }
        }
    }
    delete[] Indexes;
}

void compressibleHyperelasticity::stiffnessRhs_inhomogeneousForcing(Epetra_Vector & u, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    int node, e_gid, error;
    int n_gauss_points = Mesh->n_gauss_cells;
    double gauss_weight;
    
    int *Indexes;
    Indexes = new int [3*Mesh->el_type];
    
    Epetra_SerialDenseMatrix Ke(3*Mesh->el_type,3*Mesh->el_type);
    Epetra_SerialDenseVector Re(3*Mesh->el_type);
    
    Epetra_SerialDenseMatrix deformation_gradient(3,3);
    Epetra_SerialDenseMatrix tangent_matrix(6,6);
    Epetra_SerialDenseVector scd_piola_stress(6);
    Epetra_SerialDenseMatrix block_scd_piola_stress(9,9);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
    Epetra_SerialDenseMatrix matrix_X(3,Mesh->el_type);
    Epetra_SerialDenseMatrix xg(3,n_gauss_points);
    
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
    Epetra_SerialDenseMatrix B_times_TM(3*Mesh->el_type,6);
    Epetra_SerialDenseMatrix matrix_BG(9,3*Mesh->el_type);
    Epetra_SerialDenseMatrix BG_times_BSPS(3*Mesh->el_type,9);
    
    Epetra_SerialDenseVector fevol(3*Mesh->el_type);
    Epetra_SerialDenseVector fvol(3);
    
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
            for (int iddl=0; iddl<3; ++iddl){
                Re(3*inode+iddl) = 0.0;
                fevol(3*inode+iddl) = 0.0;
                Indexes[3*inode+iddl] = 3*node+iddl;
                for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                    for (int jddl=0; jddl<3; ++jddl){
                        Ke(3*inode+iddl,3*jnode+jddl) = 0.0;
                    }
                }
            }
        }
        xg.Multiply('N','N',1.0,matrix_X,Mesh->N_tetra,0.0);
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight = Mesh->gauss_weight_cells(gp);
            fvol = get_forcing(xg(0,gp),xg(1,gp),xg(2,gp),e_lid,gp);
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                dx_shape_functions(inode,0) = Mesh->DX_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,1) = Mesh->DY_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,2) = Mesh->DZ_N_tetra(gp+n_gauss_points*inode,e_lid);
                for (unsigned int iddl=0; iddl<3; ++iddl){
                    fevol(3*inode+iddl) += gauss_weight*pressure_load*fvol(iddl)*Mesh->N_tetra(inode,gp)*Mesh->detJac_tetra(e_lid,gp);
                }
            }
            
            error = deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
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
            
            error = Re.Multiply('T','N',-gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,scd_piola_stress,1.0);
            
            error = B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,tangent_matrix,0.0);
            error = Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
            
            error = BG_times_BSPS.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_BG,block_scd_piola_stress,0.0);
            error = Ke.Multiply('N','N',1.0,BG_times_BSPS,matrix_BG,1.0);
        }
        
        for (unsigned int i=0; i<3*Mesh->el_type; ++i){
            error = F.SumIntoGlobalValues(1, &Indexes[i], &Re(i));
            error = F.SumIntoGlobalValues(1, &Indexes[i], &fevol(i));
            for (unsigned int j=0; j<3*Mesh->el_type; ++j){
                error = K.SumIntoGlobalValues(1, &Indexes[i], 1, &Indexes[j], &Ke(i,j));
            }
        }
    }
    delete[] Indexes;
}

void compressibleHyperelasticity::rhs_NeumannBoundaryCondition(Epetra_FEVector & F){
    
    int node;
    int* Indexes;
    unsigned int e_gid;
    Indexes = new int [3*Mesh->face_type];
    
    int n_gauss_points = Mesh->n_gauss_faces;
    double gauss_weight;
    
    Epetra_SerialDenseMatrix xg(3,n_gauss_points), matrix_X(3,Mesh->face_type);
    Epetra_SerialDenseVector force(3*Mesh->face_type);
    Epetra_SerialDenseVector dead_pressure(3);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_faces; ++e_lid){
        e_gid  = Mesh->local_faces[e_lid];
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            node = Mesh->faces_nodes[Mesh->face_type*e_gid+inode];
            Indexes[3*inode+0] = 3*node+0;
            Indexes[3*inode+1] = 3*node+1;
            Indexes[3*inode+2] = 3*node+2;
            matrix_X(0,inode) = Mesh->nodes_coord[3*node+0];
            matrix_X(1,inode) = Mesh->nodes_coord[3*node+1];
            matrix_X(2,inode) = Mesh->nodes_coord[3*node+2];
            for (unsigned int iddl=0; iddl<3; ++iddl){
                force(3*inode+iddl) = 0.0;
            }
        }
        xg.Multiply('N','T',1.0,matrix_X,Mesh->N_faces,0.0);
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight  = Mesh->gauss_weight_faces(gp);
            dead_pressure = get_neumannBc(matrix_X,xg,gp);
            for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
                for (unsigned int iddl=0; iddl<3; ++iddl){
                    force(3*inode+iddl) += gauss_weight*pressure_load*dead_pressure(iddl)*Mesh->N_faces(gp,inode);
                }
            }
        }
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            for (unsigned int iddl=0; iddl<3; ++iddl){
                F.SumIntoGlobalValues(1, &Indexes[3*inode+iddl], &force(3*inode+iddl));
            }
        }
    }
    delete[] Indexes;
}
