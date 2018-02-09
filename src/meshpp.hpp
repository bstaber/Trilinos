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
    mesh(std::string & fileName_mesh, double scaling);
    mesh(Epetra_Comm & comm, std::string & fileName_mesh, double scaling);
    ~mesh();

    int read_gmsh(std::string & fileName_mesh, double scaling);
    int read_boundary_file(std::string & fileName_bc, unsigned int & number_physical_groups);
    int metis_part_mesh(int & NumProc);
    void get_local_nodes(int & MyPID);
    void get_cells_and_ghosts(int & MyPID);
    void get_cartesian_coordinate(unsigned int & e_gid, unsigned int & gp);

    void store_feinterp_faces();
    void store_feinterp_cells();

    Epetra_Comm* Comm;

    Epetra_SerialDenseMatrix N_faces, D1_N_faces, D2_N_faces;

    Epetra_SerialDenseVector local_rows, vol_cells, N_cells, detJac_cells, DX_N_cells, DY_N_cells, DZ_N_cells;

    Epetra_IntSerialDenseMatrix nodes_to_boundaries;

    std::vector<double> nodes_coord;
    std::vector<int> cells_nodes, faces_nodes;
    std::vector<int> local_nodes_without_ghosts, local_dof_without_ghosts;
    std::vector<int> local_nodes, local_dof;
    std::vector<int> local_cells, local_faces;

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
    Epetra_SerialDenseVector gauss_weight_cells, gauss_weight_faces;
    Epetra_SerialDenseVector xi_cells, eta_cells, zeta_cells;
    Epetra_SerialDenseVector xi_faces, eta_faces;
};

#endif
