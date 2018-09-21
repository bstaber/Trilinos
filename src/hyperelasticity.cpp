/*
Brian Staber (brian.staber@gmail.com)
*/

#include "hyperelasticity.hpp"
#include "fepp.hpp"

hyperelasticity::hyperelasticity(){
}

hyperelasticity::~hyperelasticity(){
}

void hyperelasticity::create_FECrsGraph(){
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

int hyperelasticity::compute_green_lagrange(Epetra_Vector & x, double & xi, double & eta, double & zeta, std::string & filename){

    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);

    Epetra_Map    CellsMap(-1,Mesh->n_local_cells,&Mesh->local_cells[0],0,*Comm);
    Epetra_Vector green_lagrange(CellsMap);

    int node, e_gid;
    double det_jac_cells;

    Epetra_SerialDenseMatrix deformation_gradient(3,3);
    Epetra_SerialDenseMatrix right_cauchy(3,3);
    Epetra_SerialDenseMatrix matrix_X(3,Mesh->el_type);
    Epetra_SerialDenseMatrix matrix_x(3,Mesh->el_type);
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix D(Mesh->el_type,3), DX(Mesh->el_type,3);
    Epetra_SerialDenseMatrix JacobianMatrix(3,3);

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
    };

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

        jacobian_matrix(matrix_X,D,JacobianMatrix);
        jacobian_det(JacobianMatrix,det_jac_cells);
        dX_shape_functions(D,JacobianMatrix,det_jac_cells,dx_shape_functions);

        deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
        right_cauchy.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);
        green_lagrange[e_lid] = 0.5*(right_cauchy(1,1)-1.0);
    }

    int NumTargetElements = 0;
    if (Comm->MyPID()==0){
        NumTargetElements = Mesh->n_cells;
    }
    Epetra_Map         MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export      ExportOnRoot(CellsMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(green_lagrange,ExportOnRoot,Insert);

    int error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
    return error;
}

void hyperelasticity::compute_center_cauchy_stress(Epetra_Vector & x, std::string & filename){

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
    double det_jac_cells, det, gauss_weight, theta;

    Epetra_SerialDenseMatrix deformation_gradient(3,3);
    Epetra_SerialDenseMatrix right_cauchy(3,3);
    Epetra_SerialDenseMatrix matrix_X(3,Mesh->el_type), matrix_x(3,Mesh->el_type);
    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3), D(Mesh->el_type,3);
    Epetra_SerialDenseMatrix JacobianMatrix(3,3);
    Epetra_SerialDenseMatrix piola_stress(3,3);
    Epetra_SerialDenseMatrix cauchy_stress(3,3);
    Epetra_SerialDenseMatrix dg_times_ps(3,3);

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
            matrix_X(0,inode) = Mesh->nodes_coord[3*node+0];
            matrix_X(1,inode) = Mesh->nodes_coord[3*node+1];
            matrix_X(2,inode) = Mesh->nodes_coord[3*node+2];
            matrix_x(0,inode) = u[OverlapMap->LID(3*node+0)] + Mesh->nodes_coord[3*node+0];
            matrix_x(1,inode) = u[OverlapMap->LID(3*node+1)] + Mesh->nodes_coord[3*node+1];
            matrix_x(2,inode) = u[OverlapMap->LID(3*node+2)] + Mesh->nodes_coord[3*node+2];
        }

        jacobian_matrix(matrix_X,D,JacobianMatrix);
        jacobian_det(JacobianMatrix,det_jac_cells);
        dX_shape_functions(D,JacobianMatrix,det_jac_cells,dx_shape_functions);

        deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
        get_material_parameters_for_recover(e_lid);
        get_stress_for_recover(deformation_gradient, det, piola_stress);
        dg_times_ps.Multiply('N','N',1.0,deformation_gradient,piola_stress,0.0);
        cauchy_stress.Multiply('N','T',1.0,dg_times_ps,deformation_gradient,0.0);

        sigma11[e_lid] = cauchy_stress(0,0);
        sigma22[e_lid] = cauchy_stress(1,1);
        sigma33[e_lid] = cauchy_stress(2,2);
        sigma12[e_lid] = cauchy_stress(0,1);
        sigma13[e_lid] = cauchy_stress(0,2);
        sigma23[e_lid] = cauchy_stress(1,2);
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

void hyperelasticity::compute_gauss_vonmises(Epetra_Vector & x, std::string & filename){

    Epetra_Vector u(*OverlapMap);
    u.Import(x, *ImportToOverlapMap, Insert);

    int node, e_gid;
    int n_gauss_points = Mesh->n_gauss_cells;
    double det_jac_cells, det, gauss_weight;

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
    Epetra_Vector vonmises(CellsMap);

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
                dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,e_lid);
                dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,e_lid);
            }

            deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);

            get_material_parameters(e_lid, gp);
            get_stress_for_recover(deformation_gradient, det, piola_stress);
            dg_times_ps.Multiply('N','N',1.0,deformation_gradient,piola_stress,0.0);
            cauchy_stress.Multiply('N','T',1.0,dg_times_ps,deformation_gradient,0.0);

            vonmises[Mesh->n_gauss_cells*e_lid+gp] = std::sqrt( (cauchy_stress(0,0)-cauchy_stress(1,1))*(cauchy_stress(0,0)-cauchy_stress(1,1)) + (cauchy_stress(1,1)-cauchy_stress(2,2))*(cauchy_stress(1,1)-cauchy_stress(2,2)) + (cauchy_stress(2,2)-cauchy_stress(0,0))*(cauchy_stress(2,2)-cauchy_stress(0,0)) + 6.0*(cauchy_stress(1,2)*cauchy_stress(1,2) + cauchy_stress(2,0)*cauchy_stress(2,0) + cauchy_stress(0,1)*cauchy_stress(0,1)) );
        }
    }

    int NumTargetElements = 0;
    if (Comm->MyPID()==0){
        NumTargetElements = Mesh->n_cells*Mesh->n_gauss_cells;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export ExportOnRoot(CellsMap,MapOnRoot);

    Epetra_MultiVector lhs_root11(MapOnRoot,true);
    lhs_root11.Export(vonmises,ExportOnRoot,Insert);
    std::string file11 = filename + "_gauss_vonmises.mtx";
    int error11 = EpetraExt::MultiVectorToMatrixMarketFile(file11.c_str(),lhs_root11,0,0,false);

}

void hyperelasticity::compute_B_matrices(Epetra_SerialDenseMatrix & F, Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B, Epetra_SerialDenseMatrix & BG){

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
