#ifndef MESHPP_HPP
#define MESHPP_HPP

//Standard includes
#include <stdio.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>
#include <iomanip>
//Epetra includes
#include "Epetra_MpiComm.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Time.h"
//ML includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"
#include "ml_include.h"
//EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"
//Teuchos includes
#include "Teuchos_RCP.hpp"
//Amesos includes
#include "Amesos.h"
//AztecOO includes
#include "AztecOO.h"
#include "Amesos_BaseSolver.h"
//Metis include
#include "metis.h"
//My includes
#include "fepp.hpp"

class mesh
{
    
public:
    mesh();
    mesh(std::string & fileName_mesh);
    mesh(Epetra_Comm & comm, std::string & fileName_mesh);
    ~mesh();
    
    int read_gmsh_tetra(std::string & fileName_mesh);
    int read_boundary_file(std::string & fileName_bc, unsigned int & number_physical_groups);
    int epetra_read_gmsh_tetra(std::string & fileName_msh, Epetra_Vector & cells_nodes_mpi, Epetra_Vector & nodes_coord_mpi);
    int metis_part_mesh(int & NumProc);
    void print_info();
    void get_local_nodes(int & MyPID);
    void get_cells_and_ghosts(int & MyPID);
    
    void store_feinterp_tri();
    void store_feinterp_tetra();
    
    Epetra_SerialDenseMatrix N_tri;
    Epetra_SerialDenseMatrix D1_N_tri;
    Epetra_SerialDenseMatrix D2_N_tri;
    Epetra_SerialDenseMatrix detJac_tri;
    
    Epetra_SerialDenseVector local_rows;
    Epetra_SerialDenseVector vol_tetra;
    Epetra_SerialDenseMatrix N_tetra;
    Epetra_SerialDenseMatrix detJac_tetra;
    Epetra_SerialDenseMatrix DX_N_tetra;
    Epetra_SerialDenseMatrix DY_N_tetra;
    Epetra_SerialDenseMatrix DZ_N_tetra;
    
    Epetra_IntSerialDenseMatrix nodes_to_boundaries;
    
    std::vector<double> nodes_coord;
    std::vector<int> cells_nodes;
    std::vector<int> faces_nodes;
    std::vector<int> local_nodes_without_ghosts;
    std::vector<int> local_dof_without_ghosts;
    std::vector<int> local_nodes;
    std::vector<int> local_dof;
    std::vector<int> local_cells;
    std::vector<int> local_faces;
    
    idx_t * epart;
    idx_t * npart;
    int * NumIndicesPerRow;
    
    int n_nodes = 0;
    int n_cells = 0;
    int n_faces = 0;
    int el_type = 0;
    int face_type = 0;
    int n_local_nodes_without_ghosts = 0;
    int n_local_nodes = 0;
    int n_local_cells = 0;
    int n_local_faces = 0;
    
    unsigned int n_gauss_faces;
    unsigned int n_gauss_cells;
    Epetra_SerialDenseVector gauss_weight_cells;
    Epetra_SerialDenseVector gauss_weight_faces;
    Epetra_SerialDenseVector xi_cells;
    Epetra_SerialDenseVector eta_cells;
    Epetra_SerialDenseVector zeta_cells;
    Epetra_SerialDenseVector xi_faces;
    Epetra_SerialDenseVector eta_faces;
    //double gauss_weight_tetra;
    //double gauss_weight_tri;
    
    Epetra_Comm* Comm;
};

#endif
