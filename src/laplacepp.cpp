/*
Brian Staber (brian.staber@gmail.com)
*/

#include "laplacepp.hpp"

laplace::laplace(){
}

laplace::laplace(mesh & mesh){
    Mesh = &mesh;
    Comm = Mesh->Comm;
    StandardMap        = new Epetra_Map(-1,Mesh->n_local_nodes_without_ghosts,&Mesh->local_nodes_without_ghosts[0],0,*Comm);
    OverlapMap         = new Epetra_Map(-1,Mesh->n_local_nodes,&Mesh->local_nodes[0],0,*Comm);
    ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);
    create_FECrsGraph();
}

laplace::laplace(mesh & mesh, Teuchos::ParameterList & Parameters){
    Mesh = &mesh;
    Comm = Mesh->Comm;
    std::string boundary_file = Teuchos::getParameter<std::string>(Parameters, "boundary_file");
    unsigned int number_physical_groups = Teuchos::getParameter<unsigned int>(Parameters, "nb_phys_groups");

    Mesh->read_boundary_file(boundary_file,number_physical_groups);

    StandardMap        = new Epetra_Map(-1,Mesh->n_local_nodes_without_ghosts,&Mesh->local_nodes_without_ghosts[0],0,*Comm);
    OverlapMap         = new Epetra_Map(-1,Mesh->n_local_nodes,&Mesh->local_nodes[0],0,*Comm);
    ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);
    create_FECrsGraph();
}

laplace::laplace(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
    std::string mesh_file = Teuchos::getParameter<std::string>(Parameters, "mesh_file");
    std::string boundary_file = Teuchos::getParameter<std::string>(Parameters, "boundary_file");
    unsigned int number_physical_groups = Teuchos::getParameter<unsigned int>(Parameters, "nb_phys_groups");

    Mesh = new mesh(comm, Parameters);
    Mesh->read_boundary_file(boundary_file,number_physical_groups);
    Comm = Mesh->Comm;

    StandardMap        = new Epetra_Map(-1,Mesh->n_local_nodes_without_ghosts,&Mesh->local_nodes_without_ghosts[0],0,*Comm);
    OverlapMap         = new Epetra_Map(-1,Mesh->n_local_nodes,&Mesh->local_nodes[0],0,*Comm);
    ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);
    create_FECrsGraph();
}

laplace::~laplace(){
    //delete Comm;
    //delete StandardMap;
    //delete OverlapMap;
    //delete ImportToOverlapMap;
    //delete FEGraph;
}

void laplace::create_FECrsGraph(){
    FEGraph = new Epetra_FECrsGraph(Copy,*StandardMap,100);
    int eglob, node;
    int *index;
    index = new int [Mesh->el_type];

    for (int eloc=0; eloc<Mesh->n_local_cells; ++eloc){
        eglob = Mesh->local_cells[eloc];
        for (int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
            index[inode] = node;
        }

        for (int i=0; i<Mesh->el_type; ++i){
            for (int j=0; j<Mesh->el_type; ++j){
                FEGraph->InsertGlobalIndices(1, &index[i], 1, &index[j]);
            }
        }

    }
    Comm->Barrier();
    FEGraph->GlobalAssemble();
    delete[] index;
}

void laplace::solve_aztec(Teuchos::ParameterList & Parameters, Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs, int * bc_indx, double * bc_val){

    assembling_OAZ(matrix, rhs, bc_indx, bc_val);
    lhs.PutScalar(0.0);

    Epetra_LinearProblem problem(&matrix, &lhs, &rhs);
    AztecOO solver(problem);

    solver.SetParameters(Parameters);

    solver.Iterate(2000,1e-6);

    if (Comm->MyPID()==0){
        std::cout << "laplace problem. \n";
        std::cout << "AZTEC_its \t AZTEC_res \t\t AZTEC_time \n";
        std::cout << solver.NumIters() << "\t\t" << solver.TrueResidual() << "\t\t" << solver.SolveTime() << "\n";
    }
}

void laplace::solve_amesos(Teuchos::ParameterList & Parameters, Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs, int * bc_indx, double * bc_val){

    bool display = Teuchos::getParameter<bool>(Parameters,"display");
    std::string solver_type = Teuchos::getParameter<std::string>(Parameters, "solver_type");

    /*if (Comm->MyPID()==0 && display){
    display_amesos_solvers();
    }*/

    assembling_OAZ(matrix, rhs, bc_indx, bc_val);
    lhs.PutScalar(0.0);

    Epetra_LinearProblem problem(&matrix, &lhs, &rhs);

    Amesos_BaseSolver* Solver;
    Amesos Factory;
    Solver = Factory.Create(solver_type, problem);
    Solver->SymbolicFactorization();
    Solver->NumericFactorization();
    Solver->Solve();
}

void laplace::assembling_OAZ(Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs, int * bc_indx, double * bc_val){

    assembling(matrix,rhs);

    int ind_bc1 = bc_indx[0];
    int ind_bc2 = bc_indx[1];
    double val_bc1 = bc_val[0];
    double val_bc2 = bc_val[1];

    Epetra_MultiVector v(*StandardMap,true);
    v.PutScalar(0.0);

    int n_BCDof = 0;
    int node;
    for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
        if(Mesh->nodes_to_boundaries(i,ind_bc1)==1 || Mesh->nodes_to_boundaries(i,ind_bc2)==1){
            n_BCDof+=1;
        }
    }

    int indbc = 0;
    int * dofOnBoundary = new int [n_BCDof];

    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes_without_ghosts[inode];
        if (Mesh->nodes_to_boundaries(inode,ind_bc1)==1){
            dofOnBoundary[indbc] = inode;
            v[0][StandardMap->LID(node)] = val_bc1;
            indbc++;
        }

        if (Mesh->nodes_to_boundaries(inode,ind_bc2)==1){
            dofOnBoundary[indbc] = inode;
            v[0][StandardMap->LID(node)] = val_bc2;
            indbc++;
        }
    }


    Epetra_MultiVector rhsDir(*StandardMap,true);
    matrix.Apply(v,rhsDir);
    rhs.Update(-1.0,rhsDir,1.0);

    for (int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes_without_ghosts[inode];
        if (Mesh->nodes_to_boundaries(inode,ind_bc1)==1){
            rhs[0][StandardMap->LID(node)] = v[0][StandardMap->LID(node)];
        }
        if (Mesh->nodes_to_boundaries(inode,ind_bc2)==1){
            rhs[0][StandardMap->LID(node)] = v[0][StandardMap->LID(node)];
        }
    }

    ML_Epetra::Apply_OAZToMatrix(dofOnBoundary,n_BCDof,matrix);
    delete[] dofOnBoundary;
}

void laplace::assembling(Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs){

    matrix.PutScalar(0.0);
    rhs.PutScalar(0.0);

    Epetra_SerialDenseMatrix B(3,Mesh->el_type);
    Epetra_SerialDenseMatrix Ke(Mesh->el_type,Mesh->el_type);

    double gaussWeight; // = 1.0/24.0;
    int n_gauss_points = Mesh->n_gauss_cells;
    unsigned int eglob;
    int *index = new int [Mesh->el_type];
    int error;

    for (unsigned int eloc=0; eloc<Mesh->n_local_cells; ++eloc){
        eglob = Mesh->local_cells[eloc];
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            index[inode] = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
            for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                Ke(inode,jnode) = 0.0;
            }
        }

        for (unsigned int gp=0; gp<Mesh->n_gauss_cells; ++gp){
            gaussWeight = Mesh->gauss_weight_cells(gp);
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                B(0,inode) = Mesh->DX_N_cells(gp+n_gauss_points*inode,eloc);
                B(1,inode) = Mesh->DY_N_cells(gp+n_gauss_points*inode,eloc);
                B(2,inode) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,eloc);
            }
            Ke.Multiply('T','N',gaussWeight*Mesh->detJac_cells(eloc,gp),B,B,1.0);
        }

        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                error = matrix.SumIntoGlobalValues(1, &index[inode], 1, &index[jnode], &Ke(inode,jnode));
            }
        }

    }

    Comm->Barrier();

    matrix.GlobalAssemble();
    matrix.FillComplete();
    rhs.GlobalAssemble();

    delete[] index;
}

void laplace::compute_local_directions(Epetra_Vector & laplace_one, Epetra_Vector & laplace_two){

    Epetra_Vector u_phi(*OverlapMap);
    u_phi.Import(laplace_one, *ImportToOverlapMap, Insert);

    Epetra_Vector u_psi(*OverlapMap);
    u_psi.Import(laplace_two, *ImportToOverlapMap, Insert);

    int node;
    unsigned int e_gid;
    double norm_phi, norm_psi;
    Epetra_SerialDenseMatrix matrix_B(3,Mesh->el_type);
    Epetra_SerialDenseVector phi_nodes(Mesh->el_type);
    Epetra_SerialDenseVector psi_nodes(Mesh->el_type);
    Epetra_SerialDenseVector grad_psi(3);
    Epetra_SerialDenseVector grad_phi(3);
    Epetra_SerialDenseMatrix matrix_X(3,Mesh->el_type);

    Epetra_SerialDenseVector e0(3);
    Epetra_SerialDenseVector e1(3);
    Epetra_SerialDenseVector e2(3);

    int n_gauss_points = Mesh->n_gauss_cells;
    laplace_direction_one.Reshape(n_gauss_points*Mesh->n_local_cells,3);
    laplace_direction_two.Reshape(n_gauss_points*Mesh->n_local_cells,3);
    laplace_direction_two_cross_one.Reshape(n_gauss_points*Mesh->n_local_cells,3);

    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            phi_nodes(inode) = u_phi[OverlapMap->LID(node)];
            psi_nodes(inode) = u_psi[OverlapMap->LID(node)];
            matrix_X(0,inode) = Mesh->nodes_coord[3*node+0];
            matrix_X(1,inode) = Mesh->nodes_coord[3*node+1];
            matrix_X(2,inode) = Mesh->nodes_coord[3*node+2];
        }

        for (unsigned int gp=0; gp<n_gauss_points; ++gp){
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                matrix_B(0,inode) = Mesh->DX_N_cells(gp+n_gauss_points*inode,e_lid);
                matrix_B(1,inode) = Mesh->DY_N_cells(gp+n_gauss_points*inode,e_lid);
                matrix_B(2,inode) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,e_lid);
            }

            grad_phi.Multiply('N','N',1.0,matrix_B,phi_nodes,0.0);
            grad_psi.Multiply('N','N',1.0,matrix_B,psi_nodes,0.0);

            norm_phi = grad_phi.Norm2();
            norm_psi = grad_psi.Norm2();

            e0(0) = grad_phi(0)/norm_phi;
            e0(1) = grad_phi(1)/norm_phi;
            e0(2) = grad_phi(2)/norm_phi;

            e2(0) = grad_psi(0)/norm_psi;
            e2(1) = grad_psi(1)/norm_psi;
            e2(2) = grad_psi(2)/norm_psi;

            e1(0) = e2(1)*e0(2) - e2(2)*e0(1);
            e1(1) = e2(2)*e0(0) - e2(0)*e0(2);
            e1(2) = e2(0)*e0(1) - e2(1)*e0(0);

            laplace_direction_one(n_gauss_points*e_lid+gp,0) = e0(0);
            laplace_direction_one(n_gauss_points*e_lid+gp,1) = e0(1);
            laplace_direction_one(n_gauss_points*e_lid+gp,2) = e0(2);

            laplace_direction_two(n_gauss_points*e_lid+gp,0) = e2(0);
            laplace_direction_two(n_gauss_points*e_lid+gp,1) = e2(1);
            laplace_direction_two(n_gauss_points*e_lid+gp,2) = e2(2);

            laplace_direction_two_cross_one(n_gauss_points*e_lid+gp,0) = e1(0);
            laplace_direction_two_cross_one(n_gauss_points*e_lid+gp,1) = e1(1);
            laplace_direction_two_cross_one(n_gauss_points*e_lid+gp,2) = e1(2);
        }
    }
}

void laplace::compute_center_local_directions(Epetra_Vector & laplace_one, Epetra_Vector & laplace_two){

    Epetra_Vector u_phi(*OverlapMap);
    u_phi.Import(laplace_one, *ImportToOverlapMap, Insert);

    Epetra_Vector u_psi(*OverlapMap);
    u_psi.Import(laplace_two, *ImportToOverlapMap, Insert);

    int node;
    unsigned int e_gid;
    double norm_phi, norm_psi;
    Epetra_SerialDenseVector phi_nodes(Mesh->el_type);
    Epetra_SerialDenseVector psi_nodes(Mesh->el_type);
    Epetra_SerialDenseVector grad_psi(3);
    Epetra_SerialDenseVector grad_phi(3);
    Epetra_SerialDenseMatrix matrix_X(3,Mesh->el_type);
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3), D(Mesh->el_type,3);
    Epetra_SerialDenseMatrix JacobianMatrix(3,3);

    Epetra_SerialDenseVector e0(3);
    Epetra_SerialDenseVector e1(3);
    Epetra_SerialDenseVector e2(3);

    int n_gauss_points = Mesh->n_gauss_cells;
    laplace_direction_one_center.Reshape(Mesh->n_local_cells,3);
    laplace_direction_two_center.Reshape(Mesh->n_local_cells,3);
    laplace_direction_two_cross_one_center.Reshape(Mesh->n_local_cells,3);

    double det_jac_cells;
    switch (Mesh->el_type){
        double xi;
        case 4:
            xi = 1.0/3.0;
            tetra4::d_shape_functions(D, xi, xi, xi);
            break;
        case 8:
            xi = 0.0;
            hexa8::d_shape_functions(D, xi, xi, xi);
            break;
        case 10:
            xi = 1.0/3.0;
            tetra10::d_shape_functions(D, xi, xi, xi);
            break;
    };
    for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        e_gid = Mesh->local_cells[e_lid];
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            phi_nodes(inode) = u_phi[OverlapMap->LID(node)];
            psi_nodes(inode) = u_psi[OverlapMap->LID(node)];
            matrix_X(0,inode) = Mesh->nodes_coord[3*node+0];
            matrix_X(1,inode) = Mesh->nodes_coord[3*node+1];
            matrix_X(2,inode) = Mesh->nodes_coord[3*node+2];
        }

        jacobian_matrix(matrix_X,D,JacobianMatrix);
        jacobian_det(JacobianMatrix,det_jac_cells);
        dX_shape_functions(D,JacobianMatrix,det_jac_cells,dx_shape_functions);

        grad_phi.Multiply('T','N',1.0,dx_shape_functions,phi_nodes,0.0);
        grad_psi.Multiply('T','N',1.0,dx_shape_functions,psi_nodes,0.0);

        norm_phi = grad_phi.Norm2();
        norm_psi = grad_psi.Norm2();

        e0(0) = grad_phi(0)/norm_phi;
        e0(1) = grad_phi(1)/norm_phi;
        e0(2) = grad_phi(2)/norm_phi;

        e2(0) = grad_psi(0)/norm_psi;
        e2(1) = grad_psi(1)/norm_psi;
        e2(2) = grad_psi(2)/norm_psi;

        e1(0) = e2(1)*e0(2) - e2(2)*e0(1);
        e1(1) = e2(2)*e0(0) - e2(0)*e0(2);
        e1(2) = e2(0)*e0(1) - e2(1)*e0(0);

        laplace_direction_one_center(e_lid,0) = e0(0);
        laplace_direction_one_center(e_lid,1) = e0(1);
        laplace_direction_one_center(e_lid,2) = e0(2);

        laplace_direction_two_center(e_lid,0) = e2(0);
        laplace_direction_two_center(e_lid,1) = e2(1);
        laplace_direction_two_center(e_lid,2) = e2(2);

        laplace_direction_two_cross_one_center(e_lid,0) = e1(0);
        laplace_direction_two_cross_one_center(e_lid,1) = e1(1);
        laplace_direction_two_cross_one_center(e_lid,2) = e1(2);
    }

}

int laplace::print_solution(Epetra_Vector & lhs, std::string fileName){
    int NumTargetElements = 0;
    if (Comm->MyPID()==0){
        NumTargetElements = Mesh->n_nodes;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export ExportOnRoot(*StandardMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(lhs,ExportOnRoot,Insert);
    int error = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(),lhs_root,0,0,false);
    return error;
}
