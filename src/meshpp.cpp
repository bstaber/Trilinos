/*
Brian Staber (brian.staber@gmail.com)
*/

#include "meshpp.hpp"
#include "error_messager.hpp"

mesh::mesh(){
}

mesh::mesh(std::string & fileName_mesh, double scaling){
    read_gmsh(fileName_mesh, scaling);
}

mesh::mesh(Epetra_Comm & comm, Teuchos::ParameterList & paramList) {
    //std::string & fileName_mesh, double scaling){
    Comm = &comm;
    int MyPID = Comm->MyPID();
    int NumProc = Comm->NumProc();

    std::string fileName_mesh = Teuchos::getParameter<std::string>(paramList.sublist("Mesh"), "mesh_file");
    double scaling = Teuchos::getParameter<double>(paramList.sublist("Mesh"), "scaling");
    MeshParameters = &(paramList.sublist("Mesh"));

    read_gmsh(fileName_mesh, scaling);

    epart = new idx_t[n_cells];
    npart = new idx_t[n_nodes];

    Comm->Barrier();
    if (Comm->NumProc()>1){
        if (MyPID==0) metis_part_mesh(NumProc);
    }
    else{
        for (unsigned int e=0; e<n_cells; ++e) epart[e] = 0;
        for (unsigned int n=0; n<n_nodes; ++n) npart[n] = 0;
    }

    Comm->Broadcast(epart,n_cells,0);
    Comm->Broadcast(npart,n_nodes,0);
    Comm->Barrier();
    get_local_nodes(MyPID);
    get_cells_and_ghosts(MyPID);

    if (n_local_faces>0 && (face_type==3 || face_type==4 || face_type==6)) store_feinterp_faces();
    store_feinterp_cells();

    Comm->Barrier();
    if (MyPID==0){
        std::cout << std::setw(5) << "MyPID" << std::setw(20) << "n_cells" << std::setw(20) << "el_type" << std::setw(20) << "face_type" << std::setw(20) << "n_nodes" << std::setw(20) << "n_faces" << std::setw(20) << "n_gauss_cells" << std::setw(20) << "n_gauss_faces" << std::setw(20) << "processors\n";
        std::cout << std::setw(5) << MyPID << std::setw(20) << n_cells << std::setw(20) << el_type << std::setw(20) << face_type << std::setw(20) << n_nodes << std::setw(20) << n_faces << std::setw(20) << n_gauss_cells << std::setw(20) << n_gauss_faces << std::setw(20) << Comm->NumProc() << "\n";
        std::cout << std::setw(5) << "MyPID" << std::setw(20) << "n_local_cells" << std::setw(20) << "n_local_nodes" << std::setw(20) << "n_local_faces" << "\n";
    }
    Comm->Barrier();
    std::cout << std::setw(5) << MyPID << std::setw(20) << n_local_cells << std::setw(20) << n_local_nodes_without_ghosts << std::setw(20) << n_local_faces << "\n";
}

mesh::~mesh(){
    delete[] epart;
    delete[] npart;
}

int mesh::read_gmsh(std::string & fileName_mesh, double scaling){
    int error = 0;
    std::ifstream meshfile;
    meshfile.open(fileName_mesh.c_str());

    if (!meshfile){
        std::cout << "*ERR* Can't open the input file\n ";
        error = 1;
        return error;
    }

    char buf[100],c;
    double xyz[3];
    unsigned int n_total_cells;
    unsigned int el_info;
    unsigned int num;
    unsigned int nbtag;
    unsigned int tag1;
    unsigned int tag2;
    unsigned int n_faces3 = 0;
    unsigned int n_quad4 = 0;
    unsigned int n_faces6 = 0;
    unsigned int n_cells4 = 0;
    unsigned int n_cells8 = 0;
    unsigned int n_cells10 = 0;
    unsigned int n_cells20 = 0;
    unsigned int n_cells27 = 0;
    unsigned int node;
    unsigned int nodes_faces3[3];
    unsigned int nodes_quad4[4];
    unsigned int nodes_faces6[6];
    unsigned int nodes_cells4[4];
    unsigned int nodes_cells8[8];
    unsigned int nodes_cells10[10];
    unsigned int nodes_cells20[20];
    unsigned int nodes_cells27[27];

    std::vector<int> tri3_nodes;
    std::vector<int> quad4_nodes, quad4_phases;
    std::vector<int> tri6_nodes;
    std::vector<int> tetra4_nodes;
    std::vector<int> hexa8_nodes;
    std::vector<int> tetra10_nodes;
    std::vector<int> hexa20_nodes;
    std::vector<int> hexa27_nodes;

    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile >> n_nodes;
    meshfile.get(c);

    nodes_coord.reserve(3*n_nodes);

    for (int i=0; i<n_nodes; ++i){
        meshfile >> num;
        meshfile >> xyz[0];
        meshfile >> xyz[1];
        meshfile >> xyz[2];
        nodes_coord[3*i+0] = xyz[0]/scaling;
        nodes_coord[3*i+1] = xyz[1]/scaling;
        nodes_coord[3*i+2] = xyz[2]/scaling;
    }

    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile >> n_total_cells;
    meshfile.getline(buf,100);

    for (int i=0; i<n_total_cells; ++i){
        meshfile >> num;
        meshfile >> el_info;
        meshfile >> nbtag;
        meshfile >> tag1;
        meshfile >> tag2;

        switch (el_info) {
            case 1:
                for (unsigned int inode=0; inode<2; ++inode){
                    meshfile >> node;
                }
                //std::cout << "2-node line not supported" << std::endl;
                //NotImplementedException();
                break;
            case 2:
                for (unsigned int inode=0; inode<3; ++inode){
                    meshfile >> nodes_faces3[inode];
                    tri3_nodes.push_back(nodes_faces3[inode]-1);
                }
                break;
            case 3:
                for (unsigned int inode=0; inode<4; ++inode){
                    meshfile >> nodes_quad4[inode];
                    quad4_nodes.push_back(nodes_quad4[inode]-1);
                }
                quad4_phases.push_back(tag1);
                break;
            case 4:
                for (unsigned int inode=0; inode<4; ++inode){
                    meshfile >> nodes_cells4[inode];
                    tetra4_nodes.push_back(nodes_cells4[inode]-1);
                }
                break;
            case 5:
                for (unsigned int inode=0; inode<8; ++inode){
                    meshfile >> nodes_cells8[inode];
                    hexa8_nodes.push_back(nodes_cells8[inode]-1);
                }
                break;
            case 9:
                for (unsigned int inode=0; inode<6; ++inode){
                    meshfile >> nodes_faces6[inode];
                    tri6_nodes.push_back(nodes_faces6[inode]-1);
                }
                break;
            case 10:
                for (unsigned int inode=0; inode<9; ++inode){
                    meshfile >>node;
                }
                //std::cout << "9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face) not supported" << std::endl;
                //NotImplementedException();
                break;
            case 11:
                for (unsigned int inode=0; inode<10; ++inode){
                    meshfile >> nodes_cells10[inode];
                    tetra10_nodes.push_back(nodes_cells10[inode]-1);
                }
                break;
            case 12:
                for (unsigned int inode=0; inode<27; ++inode){
                   meshfile >> nodes_cells27[inode];
                   hexa27_nodes.push_back(nodes_cells27[inode]-1);
                }
                break;
            case 15:
                for (unsigned int inode=0; inode<1; ++inode){
                    meshfile >> node;
                }
                //std::cout << "1-node point not supported" << std::endl;
                //NotImplementedException();
                break;
            case 16:
                for (unsigned int inode=0; inode<8; ++inode){
                    meshfile >> node;
                }
                //std::cout << "8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges) not supported" << std::endl;
                //NotImplementedException();
                break;
            case 17:
                for (unsigned int inode=0; inode<20; ++inode){
                    meshfile >> nodes_cells20[inode];
                    hexa20_nodes.push_back(nodes_cells20[inode]-1);
                }
                break;
            default:
                std::cout << "Element not supported encountered: el_info = " << el_info << "\n";
                //NotImplementedException();
                break;
        };
    }
    meshfile.close();

    n_faces3 = tri3_nodes.size()/3;
    n_quad4 = quad4_nodes.size()/4;
    n_faces6 = tri6_nodes.size()/6;
    n_cells4 = tetra4_nodes.size()/4;
    n_cells8 = hexa8_nodes.size()/8;
    n_cells10 = tetra10_nodes.size()/10;
    n_cells20 = hexa20_nodes.size()/20;
    n_cells27 = hexa27_nodes.size()/27;

    if (n_cells4==0 && n_cells10==0 && n_cells8==0 && n_cells20==0){
        std::cerr << "Your mesh is empty!\n";}
        //NotImplementedException();
    if ( (n_cells4>0 && n_cells10>0) || (n_cells4>0 && n_cells8>0) || (n_cells10>0 && n_cells8>0)){
        std::cerr << "We do not handle mixed meshes that contain tetra4's and/or tetra10's and/or hexa8's.\n";
        //NotImplementedException();
        error = 1;
        return error;
    }
    if ( (n_faces3>0 && n_faces6>0) || (n_faces3>0 && n_quad4>0) || (n_faces6>0 && n_quad4>0) ){
        std::cerr << "We do not handle mixed meshes that contain tri3's and/or tri6's and/or quad4's.\n";
        //NotImplementedException();
        error = 1;
        return error;
    }
    if (n_faces3>0){
        face_type = 3;
        n_faces = n_faces3;
        faces_nodes.reserve(tri3_nodes.size());
        faces_nodes = tri3_nodes;
        if (MeshParameters->isParameter("n_int_points_faces")) {
          int n_int_points = Teuchos::getParameter<int>(*MeshParameters, "n_int_points_faces");
          switch (n_int_points) {
            case 1:
              gauss_points_tri1(gauss_weight_faces,xi_faces,eta_faces);
              break;
            case 3:
              gauss_points_tri3(gauss_weight_faces,xi_faces,eta_faces);
              break;
            case 4:
              gauss_points_tri4(gauss_weight_faces,xi_faces,eta_faces);
              break;
            default:
              gauss_points_tri3(gauss_weight_faces,xi_faces,eta_faces);
          }
        }
        else{
          gauss_points_tri3(gauss_weight_faces,xi_faces,eta_faces);
        }
        n_gauss_faces = gauss_weight_faces.Length();
    }
    if (n_faces6>0){
        face_type = 6;
        n_faces = n_faces6;
        faces_nodes.reserve(tri6_nodes.size());
        faces_nodes = tri6_nodes;
        if (MeshParameters->isParameter("n_int_points_faces")) {
          int n_int_points = Teuchos::getParameter<int>(*MeshParameters, "n_int_points_faces");
          switch (n_int_points) {
            case 1:
              gauss_points_tri1(gauss_weight_faces,xi_faces,eta_faces);
              break;
            case 3:
              gauss_points_tri3(gauss_weight_faces,xi_faces,eta_faces);
              break;
            case 4:
              gauss_points_tri4(gauss_weight_faces,xi_faces,eta_faces);
              break;
            default:
              gauss_points_tri4(gauss_weight_faces,xi_faces,eta_faces);
          }
        }
        else{
          gauss_points_tri4(gauss_weight_faces,xi_faces,eta_faces);
        }
        n_gauss_faces = gauss_weight_faces.Length();
    }
    if (n_quad4>0){
        face_type = 4;
        n_faces = n_quad4;
        faces_nodes.reserve(quad4_nodes.size());
        faces_nodes = quad4_nodes;
        faces_phases = quad4_phases;

        if (MeshParameters->isParameter("n_int_points_faces")) {
          int n_int_points = Teuchos::getParameter<int>(*MeshParameters, "n_int_points_faces");
          switch (n_int_points) {
            case 1:
              gauss_points_quad1(gauss_weight_faces,xi_faces,eta_faces);
              break;
            case 4:
              gauss_points_quad4(gauss_weight_faces,xi_faces,eta_faces);
              break;
            case 9:
              gauss_points_quad9(gauss_weight_faces,xi_faces,eta_faces);
              break;
            default:
              gauss_points_quad4(gauss_weight_faces,xi_faces,eta_faces);
          }
        }
        else{
          gauss_points_quad4(gauss_weight_faces,xi_faces,eta_faces);
        }
        n_gauss_faces = gauss_weight_faces.Length();
    }
    if (n_cells4>0){
        n_cells = n_cells4;
        el_type = 4;
        cells_nodes.reserve(tetra4_nodes.size());
        cells_nodes = tetra4_nodes;

        if (MeshParameters->isParameter("n_int_points_cells")) {
          int n_int_points = Teuchos::getParameter<int>(*MeshParameters, "n_int_points_cells");
          switch (n_int_points) {
            case 1:
              gauss_points_tetra1(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 4:
              gauss_points_tetra4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 5:
              gauss_points_tetra5(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 11:
              gauss_points_tetra11(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            default:
              gauss_points_tetra4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
          }
        }
        else{
          gauss_points_tetra4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
        }
        n_gauss_cells = gauss_weight_cells.Length();
    }
    if (n_cells8>0){
        n_cells = n_cells8;
        el_type = 8;
        cells_nodes.reserve(hexa8_nodes.size());
        cells_nodes = hexa8_nodes;

        if (MeshParameters->isParameter("n_int_points_cells")) {
          int n_int_points = Teuchos::getParameter<int>(*MeshParameters, "n_int_points_cells");
          switch (n_int_points) {
            case 4:
              gauss_points_hexa4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 8:
              gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 27:
              gauss_points_hexa27(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            default:
              gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
          }
        }
        else{
          gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
        }
        n_gauss_cells = gauss_weight_cells.Length();
    }
    if (n_cells10>0){
        n_cells = n_cells10;
        el_type= 10;
        cells_nodes.reserve(tetra10_nodes.size());
        cells_nodes = tetra10_nodes;

        if (MeshParameters->isParameter("n_int_points_cells")) {
          int n_int_points = Teuchos::getParameter<int>(*MeshParameters, "n_int_points_cells");
          switch (n_int_points) {
            case 1:
              gauss_points_tetra1(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 4:
              gauss_points_tetra4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 5:
              gauss_points_tetra5(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 11:
              gauss_points_tetra11(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            default:
              gauss_points_tetra4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
          }
        }
        else{
          gauss_points_tetra11(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
        }
        n_gauss_cells = gauss_weight_cells.Length();
    }
    if (n_cells20>0){
        n_cells = n_cells20;
        el_type = 20;
        cells_nodes.reserve(hexa20_nodes.size());
        cells_nodes = hexa20_nodes;

        if (MeshParameters->isParameter("n_int_points_cells")) {
          int n_int_points = Teuchos::getParameter<int>(*MeshParameters, "n_int_points_cells");
          switch (n_int_points) {
            case 4:
              gauss_points_hexa4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 8:
              gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 27:
              gauss_points_hexa27(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            default:
              gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
          }
        }
        else{
          gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
        }
        n_gauss_cells = gauss_weight_cells.Length();
        //std::cout << "HEXA20 SHAPE FUNCTIONS NOT IMPLEMENTED YET" << std::endl;
        //NotImplementedException();
    }
    if (n_cells27>0){
        n_cells = n_cells27;
        el_type = 27;
        cells_nodes.reserve(hexa27_nodes.size());
        cells_nodes = hexa27_nodes;

        if (MeshParameters->isParameter("n_int_points_cells")) {
          int n_int_points = Teuchos::getParameter<int>(*MeshParameters, "n_int_points_cells");
          switch (n_int_points) {
            case 4:
              gauss_points_hexa4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 8:
              gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            case 27:
              gauss_points_hexa27(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
              break;
            default:
              gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
          }
        }
        else{
          gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
        }
        n_gauss_cells = gauss_weight_cells.Length();
        std::cout << "HEXA27 SHAPE FUNCTIONS NOT IMPLEMENTED YET" << std::endl;
        //NotImplementedException();
    }
    Comm->Barrier();
    return error;
}

int mesh::read_boundary_file(std::string & fileName_bc, unsigned int & number_physical_groups){

    std::ifstream file_bc;
    file_bc.open(fileName_bc.c_str());
    if (!file_bc){
        std::cout << "*ERR* Can't open the input file called " << fileName_bc.c_str()  << "\n";
        return 1;
    }

    Epetra_IntSerialDenseVector input(number_physical_groups*n_nodes);
    for (unsigned int i=0; i<number_physical_groups*n_nodes; ++i){
        file_bc >> input(i);
    }

    nodes_to_boundaries.Reshape(n_local_nodes_without_ghosts,number_physical_groups);
    int node;
    for (unsigned int inode=0; inode<n_local_nodes_without_ghosts; ++inode){
        node = local_nodes[inode];
        for (unsigned pg=0; pg<number_physical_groups; ++pg){
        nodes_to_boundaries(inode,pg) = input(number_physical_groups*node+pg);
        }
    }
    file_bc.close();
    return 0;

}

int mesh::metis_part_mesh(int & NumProc){
    int check_PartMeshDual;

    idx_t ne = n_cells;
    idx_t nn = n_nodes;
    idx_t *eptr, *eind;
    eptr = new idx_t[n_cells+1];
    eind = new idx_t[el_type*n_cells];

    for (unsigned int i=0; i<n_cells; ++i){
        for (unsigned int inode=0; inode<el_type; ++inode){
            eind[el_type*i+inode] = cells_nodes[el_type*i+inode];
        }
        eptr[i] = el_type*i;
    }
    eptr[n_cells] = el_type*n_cells;

    idx_t *vwgt=NULL;
    idx_t *vsize=NULL;
    idx_t common;
    switch (el_type){
        case 4:
            common = 3;
            break;
        case 8:
            common = 4;
            break;
        case 10:
            common = 6;
            break;
        case 20:
            common = 8;
            break;
        default:
            std::cout << "el_type " << el_type << "not supported yet" << std::endl;
    };
    idx_t nparts = NumProc;
    real_t *tpwgts=NULL;
    idx_t *options=NULL;
    idx_t objval;
    check_PartMeshDual = METIS_PartMeshDual(&ne, &nn, eptr, eind, vwgt, vsize, &common, &nparts, tpwgts, options, &objval, epart, npart);

    if (check_PartMeshDual==0){
        std::cout << "*ERR* An error occured with METIS.\n";
    }
    delete [] eptr;
    delete [] eind;
    delete [] vwgt;
    delete [] vsize;
    delete [] tpwgts;
    delete [] options;

    return check_PartMeshDual;
}

void mesh::get_local_nodes(int & MyPID){
    for (unsigned int j=0; j<n_nodes; ++j){
        if (npart[j]==MyPID){
            local_nodes_without_ghosts.push_back(j);
            local_dof_without_ghosts.push_back(3*j+0);
            local_dof_without_ghosts.push_back(3*j+1);
            local_dof_without_ghosts.push_back(3*j+2);
        }
    }
    n_local_nodes_without_ghosts = local_nodes_without_ghosts.size();
}

void mesh::get_cells_and_ghosts(int & MyPID){

    Epetra_IntSerialDenseVector mynpart(Copy,npart,n_nodes);
    local_nodes = local_nodes_without_ghosts;
    local_dof = local_dof_without_ghosts;
    int node;

    for (unsigned int i=0; i<n_cells; ++i){
        if (epart[i]==MyPID){
            local_cells.push_back(i);
            for (unsigned inode=0; inode<el_type; ++inode){
                node = cells_nodes[el_type*i+inode];
                if (mynpart[node]!=MyPID){
                    mynpart[node]=MyPID;
                    local_nodes.push_back(node);
                    local_dof.push_back(3*node+0);
                    local_dof.push_back(3*node+1);
                    local_dof.push_back(3*node+2);
                }
            }
        }
    }

    n_local_nodes = local_nodes.size();
    n_local_cells = local_cells.size();

    if (n_faces>0){
    int nodes[face_type];
    for (unsigned int i=0; i<n_faces; ++i){
        for (unsigned int inode=0; inode<face_type; ++inode){
            nodes[inode] = faces_nodes[face_type*i+inode];
        }
        switch (face_type){
            case 3:
                if (mynpart[nodes[0]]==MyPID && mynpart[nodes[1]]==MyPID && mynpart[nodes[2]]==MyPID){
                    local_faces.push_back(i);
                    //break;
                }
                break;
            case 4:
                if (mynpart[nodes[0]]==MyPID && mynpart[nodes[1]]==MyPID && mynpart[nodes[2]]==MyPID && mynpart[nodes[3]]==MyPID){
                    local_faces.push_back(i);
                    //break;
                }
                break;
            case 6:
                if (mynpart[nodes[0]]==MyPID && mynpart[nodes[1]]==MyPID && mynpart[nodes[2]]==MyPID && mynpart[nodes[3]]==MyPID && mynpart[nodes[4]]==MyPID && mynpart[nodes[5]]==MyPID){
                    local_faces.push_back(i);
                    //break;
                }
                break;
            default:
                std::cout << "Non supported face type encountered" << std::endl;
                //NotImplementedException();
        }
    }
    n_local_faces = local_faces.size();
    }
    int check_n_faces;
    Comm->SumAll(&n_local_faces,&check_n_faces,1);
    if (check_n_faces!=n_faces){
        if (Comm->MyPID()==0){
            std::cout << "\n";
            std::cout << std::setw(15) << "****************************************************\n";
            std::cout << std::setw(15) << "ALERT: SOME FACES MAY BELONG TO MULTIPLE PROCESSOR!!\n";
            std::cout << std::setw(15) << "SOLUTION: LOWER THE NUMBER OF PROCESSORS!!\n";
            std::cout << std::setw(15) << "****************************************************\n";
        }
    }
}

Epetra_SerialDenseVector mesh::get_cartesian_coordinate(unsigned int & e_gid, unsigned int & gp){

  int node;
  int n_local_cells = n_local_cells;
  int n_gauss_cells = n_gauss_cells;
  double xi, eta, zeta;

  Epetra_SerialDenseVector vector_x(3);
  Epetra_SerialDenseVector shape_functions(el_type);
  Epetra_SerialDenseMatrix matrix_X(3,el_type);

  for (unsigned inode=0; inode<el_type; ++inode){
      node = cells_nodes[el_type*e_gid+inode];
      matrix_X(0,inode) = nodes_coord[3*node+0];
      matrix_X(1,inode) = nodes_coord[3*node+1];
      matrix_X(2,inode) = nodes_coord[3*node+2];
  }
  xi = xi_cells(gp); eta = eta_cells(gp); zeta = zeta_cells(gp);
  switch (el_type){
      case 4:
          tetra4::shape_functions(shape_functions, xi, eta, zeta);
          break;
      case 8:
          hexa8::shape_functions(shape_functions, xi, eta, zeta);
          break;
      case 10:
          tetra10::shape_functions(shape_functions, xi, eta, zeta);
          break;
      case 20:
          hexa20::shape_functions(shape_functions, xi, eta, zeta);
          break;
      default:
          std::cout << "Non supported element type encountered" << std::endl;
          //NotImplementedException();
  }
  vector_x.Multiply('N','N',1.0,matrix_X,shape_functions,0.0);
  return vector_x;
}

void mesh::store_feinterp_faces(){

    int node, eglob;
    Epetra_SerialDenseVector N(face_type);
    Epetra_SerialDenseMatrix D(face_type,2);

    N_faces.Reshape(n_gauss_faces,face_type);
    D1_N_faces.Reshape(n_gauss_faces,face_type);
    D2_N_faces.Reshape(n_gauss_faces,face_type);

    switch (face_type){
        case 3:
            for (unsigned int gp=0; gp<n_gauss_faces; ++gp){
                tri3::shape_functions(N, xi_faces[gp], eta_faces[gp]);
                for (int inode=0; inode<face_type; ++inode){
                    N_faces(gp,inode) = N(inode);
                }
            }
            break;
        case 4:
            for (unsigned int gp=0; gp<n_gauss_faces; ++gp){
                quad4::shape_functions(N, xi_faces[gp], eta_faces[gp]);
                for (int inode=0; inode<face_type; ++inode){
                    N_faces(gp,inode) = N(inode);
                }
            }
            break;
        case 6:
            for (unsigned int gp=0; gp<n_gauss_faces; ++gp){
                tri6::shape_functions(N, xi_faces[gp], eta_faces[gp]);
                for (int inode=0; inode<face_type; ++inode){
                    N_faces(gp,inode) = N(inode);
                }
            }
            break;
        default:
            std::cout << "Non supported face type encountered" << std::endl;
            //NotImplementedException();
    };

    for (unsigned int gp=0; gp<n_gauss_faces; ++gp){
        switch (face_type){
            case 3:
                tri3::d_shape_functions(D, xi_faces[gp], eta_faces[gp]);
                break;
            case 4:
                quad4::d_shape_functions(D, xi_faces[gp], eta_faces[gp]);
                break;
            case 6:
                tri6::d_shape_functions(D, xi_faces[gp], eta_faces[gp]);
                break;
            default:
                std::cout << "Non supported element type encountered" << std::endl;
        };
        for (unsigned int inode=0; inode<face_type; ++inode){
            D1_N_faces(gp,inode) = D(inode,0);
            D2_N_faces(gp,inode) = D(inode,1);
        }
    }
}

void mesh::store_feinterp_cells(){

    int node, eglob;
    double alpha, beta;
    Epetra_SerialDenseVector N(el_type);
    Epetra_SerialDenseMatrix JacobianMatrix(3,3), InverseJacobianMatrix(3,3), X(3,el_type), D(el_type,3), DX(el_type,3);

    local_rows.Resize(3*el_type*n_local_cells);
    vol_cells.Resize(n_local_cells);
    N_cells.Reshape(el_type,n_gauss_cells);
    detJac_cells.Reshape(n_local_cells,n_gauss_cells);
    DX_N_cells.Reshape(n_gauss_cells*el_type,n_local_cells);
    DY_N_cells.Reshape(n_gauss_cells*el_type,n_local_cells);
    DZ_N_cells.Reshape(n_gauss_cells*el_type,n_local_cells);

    switch (el_type){
        case 4:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                tetra4::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_cells(inode,gp) = N(inode);
                }
            }
            break;
        case 8:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                hexa8::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_cells(inode,gp) = N(inode);
                }
            }
            break;
        case 10:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                tetra10::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_cells(inode,gp) = N(inode);
                }
            }
            break;
        case 20:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                hexa20::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                  N_cells(inode,gp) = N(inode);
                }
            }
            break;
        default:
            std::cout << "Non supported element type encountered" << std::endl;
            //NotImplementedException();
    };
    //double volume = 0.0;
    for (unsigned int eloc=0; eloc<n_local_cells; ++eloc){
        eglob = local_cells[eloc];
        for (unsigned int inode=0; inode<el_type; ++inode){
            node = cells_nodes[el_type*eglob+inode];
            X(0,inode) = nodes_coord[3*node];
            X(1,inode) = nodes_coord[3*node+1];
            X(2,inode) = nodes_coord[3*node+2];
            local_rows(3*el_type*eloc+3*inode) = 3*node;
            local_rows(3*el_type*eloc+3*inode+1) = 3*node+1;
            local_rows(3*el_type*eloc+3*inode+2) = 3*node+2;
        }

        vol_cells(eloc) = 0.0;
        for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
            switch (el_type){
                case 4:
                    tetra4::d_shape_functions(D, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                    break;
                case 8:
                    hexa8::d_shape_functions(D, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                    break;
                case 10:
                    tetra10::d_shape_functions(D, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                    break;
                case 20:
                    hexa20::d_shape_functions(D, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                    break;
                default:
                    std::cout << "Non supported element type encountered" << std::endl;
                    //NotImplementedException();
            }
            jacobian_matrix(X,D,JacobianMatrix);
            jacobian_det(JacobianMatrix,detJac_cells(eloc,gp));
            dX_shape_functions(D,JacobianMatrix,detJac_cells(eloc,gp),DX);
            vol_cells(eloc) += gauss_weight_cells(gp)*detJac_cells(eloc,gp);
            for (int inode=0; inode<el_type; ++inode){
                DX_N_cells(gp+n_gauss_cells*inode,eloc) = DX(inode,0);
                DY_N_cells(gp+n_gauss_cells*inode,eloc) = DX(inode,1);
                DZ_N_cells(gp+n_gauss_cells*inode,eloc) = DX(inode,2);
            }
        }
        //volume += vol_cells(eloc);
    }
    //std::cout << "VOLUME = " << volume << "\n";
}

void mesh::update_store_feinterp_cells(Epetra_Vector & u, Epetra_Map & OverlapMap){

    int node, eglob;
    double alpha, beta;
    Epetra_SerialDenseVector N(el_type);
    Epetra_SerialDenseMatrix JacobianMatrix(3,3), InverseJacobianMatrix(3,3), X(3,el_type), D(el_type,3), DX(el_type,3);

    switch (el_type){
        case 4:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                tetra4::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_cells(inode,gp) = N(inode);
                }
            }
            break;
        case 8:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                hexa8::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_cells(inode,gp) = N(inode);
                }
            }
            break;
        case 10:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                tetra10::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_cells(inode,gp) = N(inode);
                }
            }
            break;
        case 20:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                hexa20::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_cells(inode,gp) = N(inode);
                }
            }
            break;
        default:
            std::cout << "Non supported element type encountered" << std::endl;
            //NotImplementedException();
    };
    //double volume = 0.0;
    for (unsigned int eloc=0; eloc<n_local_cells; ++eloc){
        eglob = local_cells[eloc];
        for (unsigned int inode=0; inode<el_type; ++inode){
            node = cells_nodes[el_type*eglob+inode];
            X(0,inode) = nodes_coord[3*node+0] + u[OverlapMap.LID(3*node+0)];
            X(1,inode) = nodes_coord[3*node+1] + u[OverlapMap.LID(3*node+1)];
            X(2,inode) = nodes_coord[3*node+2] + u[OverlapMap.LID(3*node+2)];
            local_rows(3*el_type*eloc+3*inode+0) = 3*node+0;
            local_rows(3*el_type*eloc+3*inode+1) = 3*node+1;
            local_rows(3*el_type*eloc+3*inode+2) = 3*node+2;
        }

        vol_cells(eloc) = 0.0;
        for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
            switch (el_type){
                case 4:
                    tetra4::d_shape_functions(D, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                    break;
                case 8:
                    hexa8::d_shape_functions(D, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                    break;
                case 10:
                    tetra10::d_shape_functions(D, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                    break;
                case 20:
                    hexa20::d_shape_functions(D, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                    break;
                default:
                    std::cout << "Non supported element type encountered" << std::endl;
                    //NotImplementedException();
            }
            jacobian_matrix(X,D,JacobianMatrix);
            jacobian_det(JacobianMatrix,detJac_cells(eloc,gp));
            dX_shape_functions(D,JacobianMatrix,detJac_cells(eloc,gp),DX);
            vol_cells(eloc) += gauss_weight_cells(gp)*detJac_cells(eloc,gp);
            for (int inode=0; inode<el_type; ++inode){
                DX_N_cells(gp+n_gauss_cells*inode,eloc) = DX(inode,0);
                DY_N_cells(gp+n_gauss_cells*inode,eloc) = DX(inode,1);
                DZ_N_cells(gp+n_gauss_cells*inode,eloc) = DX(inode,2);
            }
        }
        //volume += vol_cells(eloc);
    }
  //std::cout << "VOLUME = " << volume << std::endl;
}
