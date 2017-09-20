#include "linearelasticity_setup_pp.hpp"
#include "fepp.hpp"

LinearizedElasticity::LinearizedElasticity(){
    //dead_pressure.Resize(3);
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

void LinearizedElasticity::assemble_dirichlet(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){

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

void LinearizedElasticity::assemble_dirichlet_dead_neumann(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
    
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

void LinearizedElasticity::material_stiffness_and_rhs_dirichlet(Epetra_FECrsMatrix & K, Epetra_FEVector & F){

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
    
    Epetra_SerialDenseVector dead_pressure_load(3);
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
                    force(3*inode+iddl) += gauss_weight*dead_pressure_load(iddl)*Mesh->N_tri(gp,inode)*Mesh->detJac_tri(e_lid,gp);
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
