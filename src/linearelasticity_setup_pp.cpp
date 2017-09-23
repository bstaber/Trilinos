#include "linearelasticity_setup_pp.hpp"
#include "fepp.hpp"

LinearizedElasticity::LinearizedElasticity(){
    dead_pressure.Resize(3);
}

LinearizedElasticity::~LinearizedElasticity(){
}

void LinearizedElasticity::create_FECrsGraph(){
    
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

void LinearizedElasticity::assemble_dirichlet(Epetra_FECrsMatrix & K){

    int error;
    
    K.PutScalar(0.0);
    
    material_stiffness_and_rhs_dirichlet(K);
    
    Comm->Barrier();
    
    error=K.GlobalAssemble();
    error=K.FillComplete();
}

void LinearizedElasticity::assemble_dirichlet_dead_neumann(Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    
    F.PutScalar(0.0);
    K.PutScalar(0.0);
        
    material_stiffness_and_rhs_dirichlet(K);
    force_dead_pressure(F);
    
    Comm->Barrier();
    
    K.GlobalAssemble();
    K.FillComplete();
    F.GlobalAssemble();
}

void LinearizedElasticity::material_stiffness_and_rhs_dirichlet(Epetra_FECrsMatrix & K){

    int node, e_gid, error;
    int n_gauss_points = Mesh->n_gauss_cells;
    double gauss_weight;

    int *Indices_tetra;
    Indices_tetra = new int [3*Mesh->el_type];
    
    Epetra_SerialDenseMatrix Ke(3*Mesh->el_type,3*Mesh->el_type);
    
    Epetra_SerialDenseMatrix tangent_matrix(6,6);
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
    Epetra_SerialDenseMatrix B_times_TM(3*Mesh->el_type,6);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            for (int iddl=0; iddl<3; ++iddl){
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
            
            compute_B_matrices(dx_shape_functions,matrix_B);
            get_elasticity_tensor(e_lid, gp, tangent_matrix);
            
            error=B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_tetra(e_lid,gp),matrix_B,tangent_matrix,0.0);
            error=Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
        }
        
        for (unsigned int i=0; i<3*Mesh->el_type; ++i){
            for (unsigned int j=0; j<3*Mesh->el_type; ++j){
                error=K.SumIntoGlobalValues(1, &Indices_tetra[i], 1, &Indices_tetra[j], &Ke(i,j));
            }
        }
    }
    
}

void LinearizedElasticity::force_dead_pressure(Epetra_FEVector & F){

    int node;
    int* Indices_tri;
    unsigned int e_gid;
    Indices_tri = new int [3*Mesh->face_type];
    
    int n_gauss_points = Mesh->n_gauss_faces;
    double gauss_weight;
    
    Epetra_SerialDenseVector force(3*Mesh->face_type);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_faces; ++e_lid){
        e_gid  = Mesh->local_faces[e_lid];
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            node = Mesh->faces_nodes[Mesh->face_type*e_gid+inode];
            Indices_tri[3*inode]   = 3*node;
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
                    force(3*inode+iddl) += gauss_weight*dead_pressure(iddl)*Mesh->N_tri(gp,inode)*Mesh->detJac_tri(e_lid,gp);
                }
            }
        }
        
        for (unsigned int inode=0; inode<Mesh->face_type; ++inode){
            for (unsigned int iddl=0; iddl<3; ++iddl){
                F.SumIntoGlobalValues(1, &Indices_tri[3*inode+iddl], &force(3*inode+iddl));
            }
        }
        
    }
}

void LinearizedElasticity::compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B){
    double factor = 1.0/std::sqrt(2.0);
    for (unsigned inode=0; inode<Mesh->el_type; ++inode){
        B(0,3*inode) = dx_shape_functions(inode,0);
        B(0,3*inode+1) = 0.0;
        B(0,3*inode+2) = 0.0;
        
        B(1,3*inode) = 0.0;
        B(1,3*inode+1) = dx_shape_functions(inode,1);
        B(1,3*inode+2) = 0.0;
        
        B(2,3*inode) = 0.0;
        B(2,3*inode+1) = 0.0;
        B(2,3*inode+2) = dx_shape_functions(inode,2);
        
        B(3,3*inode) = 0.0;
        B(3,3*inode+1) = factor*dx_shape_functions(inode,2);
        B(3,3*inode+2) = factor*dx_shape_functions(inode,1);
        
        B(4,3*inode) = factor*dx_shape_functions(inode,2);
        B(4,3*inode+1) = 0.0;
        B(4,3*inode+2) = factor*dx_shape_functions(inode,0);
        
        B(5,3*inode) = factor*dx_shape_functions(inode,1);
        B(5,3*inode+1) = factor*dx_shape_functions(inode,0);
        B(5,3*inode+2) = 0.0;
    }
}

void LinearizedElasticity::compute_mean_cauchy_stress(Epetra_Vector & x, std::string & filename){
    
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
    
    Epetra_SerialDenseVector epsilon(6);
    Epetra_SerialDenseVector cauchy_stress(6);
    Epetra_SerialDenseVector B_times_epsilon(6);
    Epetra_SerialDenseVector vector_u(3*Mesh->el_type);
    Epetra_SerialDenseMatrix tangent_matrix(6,6);
    Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            vector_u(3*inode+0) = u[OverlapMap->LID(3*node+0)];
            vector_u(3*inode+1) = u[OverlapMap->LID(3*node+1)];
            vector_u(3*inode+2) = u[OverlapMap->LID(3*node+2)];
        }
        
        theta = 0.0;
        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            gauss_weight = Mesh->gauss_weight_cells(gp);
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                dx_shape_functions(inode,0) = Mesh->DX_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,1) = Mesh->DY_N_tetra(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,2) = Mesh->DZ_N_tetra(gp+n_gauss_points*inode,e_lid);
            }
            
            compute_B_matrices(dx_shape_functions,matrix_B);
            epsilon.Multiply('N','N',1.0,matrix_B,vector_u,0.0);
            get_elasticity_tensor(e_lid, gp, tangent_matrix);
            B_times_epsilon.Multiply('N','N',1.0,matrix_B,epsilon,0.0);
            cauchy_stress.Multiply('N','N',1.0,tangent_matrix,B_times_epsilon,0.0);
            
            std::cout << cauchy_stress;
            
            sigma11[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(0);
            sigma22[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(1);
            sigma33[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(2);
            sigma12[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(5);
            sigma13[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(4);
            sigma23[e_lid] += gauss_weight*Mesh->detJac_tetra(e_lid,gp)*cauchy_stress(3);
            
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

