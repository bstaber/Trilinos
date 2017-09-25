#include "meshpp.hpp"

mesh::mesh(){
}

mesh::mesh(std::string & fileName_mesh){
    read_gmsh_tetra(fileName_mesh);
}

mesh::mesh(Epetra_Comm & comm, std::string & fileName_mesh){
    Comm = &comm;
    int MyPID = Comm->MyPID();
    int NumProc = Comm->NumProc();
    
    read_gmsh_tetra(fileName_mesh);
    
    epart = new idx_t[n_cells];
    npart = new idx_t[n_nodes];
    
    Comm->Barrier();
    if (Comm->NumProc()>1){
        if (MyPID==0){
            metis_part_mesh(NumProc);
        }
    }
    else{
        for (unsigned int e=0; e<n_cells; ++e){
            epart[e] = 0;
        }
        for (unsigned int n=0; n<n_nodes; ++n){
            npart[n] = 0;
        }
    }
    
    Comm->Broadcast(epart,n_cells,0);
    Comm->Broadcast(npart,n_nodes,0);
    Comm->Barrier();
    get_local_nodes(MyPID);
    get_cells_and_ghosts(MyPID);
    
    if (MyPID==0){
        std::cout << std::setw(5) << "MyPID" << std::setw(20) << "n_cells" << std::setw(20) << "el_type" << std::setw(20) << "face_type" << std::setw(20) << "n_nodes" << std::setw(20) << "n_faces" << std::setw(20) << "processors\n";
        std::cout << std::setw(5) << MyPID << std::setw(20) << n_cells << std::setw(20) << el_type << std::setw(20) << face_type << std::setw(20) << n_nodes << std::setw(20) << n_faces << std::setw(20) << Comm->NumProc() << "\n";
        std::cout << std::setw(5) << "MyPID" << std::setw(20) << "n_local_cells" << std::setw(20) << "n_local_nodes" << std::setw(20) << "n_local_faces" << "\n";
    }
    Comm->Barrier();
    std::cerr << std::setw(5) << MyPID << std::setw(20) << n_local_cells << std::setw(20) << n_local_nodes_without_ghosts << std::setw(20) << n_local_faces << "\n";
    
    if (n_local_faces>0 && (face_type==3 || face_type==4 || face_type==6)){
        store_feinterp_tri();
    }
    store_feinterp_tetra();
}

mesh::~mesh(){
    //delete[] NumIndicesPerRow;
    delete[] epart;
    delete[] npart;
    //delete Comm;
}

int mesh::read_gmsh_tetra(std::string & fileName_mesh){
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
    unsigned int n_tri3 = 0;
    unsigned int n_quad4 = 0;
    unsigned int n_tri6 = 0;
    unsigned int n_tetra4 = 0;
    unsigned int n_hexa8 = 0;
    unsigned int n_tetra10 = 0;
    unsigned int node;
    unsigned int nodes_tri3[3];
    unsigned int nodes_quad4[4];
    unsigned int nodes_tri6[6];
    unsigned int nodes_tetra4[4];
    unsigned int nodes_hexa8[8];
    unsigned int nodes_tetra10[10];
    
    std::vector<int> tri3_nodes;
    std::vector<int> quad4_nodes;
    std::vector<int> tri6_nodes;
    std::vector<int> tetra4_nodes;
    std::vector<int> hexa8_nodes;
    std::vector<int> tetra10_nodes;
    
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
        nodes_coord[3*i] = xyz[0]/1000.0;
        nodes_coord[3*i+1] = xyz[1]/1000.0;
        nodes_coord[3*i+2] = xyz[2]/1000.0;
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
                break;
            case 2:
                for (unsigned int inode=0; inode<3; ++inode){
                    meshfile >> nodes_tri3[inode];
                    tri3_nodes.push_back(nodes_tri3[inode]-1);
                }
                break;
            case 3:
                for (unsigned int inode=0; inode<4; ++inode){
                    meshfile >> nodes_quad4[inode];
                    if (tag1==92 || tag1==93){
                        quad4_nodes.push_back(nodes_quad4[inode]-1);
                    }
                }
                break;
            case 4:
                for (unsigned int inode=0; inode<4; ++inode){
                    meshfile >> nodes_tetra4[inode];
                    tetra4_nodes.push_back(nodes_tetra4[inode]-1);
                }
                break;
            case 5:
                for (unsigned int inode=0; inode<8; ++inode){
                    meshfile >> nodes_hexa8[inode];
                    hexa8_nodes.push_back(nodes_hexa8[inode]-1);
                }
                break;
            case 9:
                for (unsigned int inode=0; inode<6; ++inode){
                    meshfile >> nodes_tri6[inode];
                    if (tag1==1){
                        tri6_nodes.push_back(nodes_tri6[inode]-1);
                    }
                }
                break;
            case 11:
                for (unsigned int inode=0; inode<10; ++inode){
                    meshfile >> nodes_tetra10[inode];
                    tetra10_nodes.push_back(nodes_tetra10[inode]-1);
                }
                break;
            case 15:
                for (unsigned int inode=0; inode<1; ++inode){
                    meshfile >> node;
                }
                break;
            default:
                std::cout << "Element not supported encountered: el_info = " << el_info << "\n";
                break;
        };
    }
    meshfile.close();
    
    n_tri3 = tri3_nodes.size()/3;
    n_quad4 = quad4_nodes.size()/4;
    n_tri6 = tri6_nodes.size()/6;
    n_tetra4 = tetra4_nodes.size()/4;
    n_hexa8 = hexa8_nodes.size()/8;
    n_tetra10 = tetra10_nodes.size()/10;
    
    if (n_tetra4==0 && n_tetra10==0 && n_hexa8==0){
        std::cerr << "Your mesh is empty!\n";}
    if ( (n_tetra4>0 && n_tetra10>0) || (n_tetra4>0 && n_hexa8>0) || (n_tetra10>0 && n_hexa8>0) ){
        std::cerr << "We do not handle mixed meshes that contain tetra4's and/or tetra10's and/or hexa8's.\n";
        error = 1;
        return error;
    }
    if ( (n_tri3>0 && n_tri6>0) || (n_tri3>0 && n_quad4>0) || (n_tri6>0 && n_quad4>0) ){
        std::cerr << "We do not handle mixed meshes that contain tri3's and/or tri6's and/or quad4's.\n";
        error = 1;
        return error;
    }
    if (n_tri3>0){
        face_type = 3;
        n_faces = n_tri3;
        faces_nodes.reserve(tri3_nodes.size());
        faces_nodes = tri3_nodes;
        
        gauss_points_tri3(gauss_weight_faces,xi_faces,eta_faces);
        n_gauss_faces = gauss_weight_faces.Length();
    }
    if (n_tri6>0){
        face_type = 6;
        n_faces = n_tri6;
        faces_nodes.reserve(tri6_nodes.size());
        faces_nodes = tri6_nodes;
        
        gauss_points_tri4(gauss_weight_faces,xi_faces,eta_faces);
        n_gauss_faces = gauss_weight_faces.Length();
    }
    if (n_quad4>0){
        face_type = 4;
        n_faces = n_quad4;
        faces_nodes.reserve(quad4_nodes.size());
        faces_nodes = quad4_nodes;
        
        gauss_points_quad4(gauss_weight_faces,xi_faces,eta_faces);
        n_gauss_faces = gauss_weight_faces.Length();
    }
    if (n_tetra4>0){
        n_cells = n_tetra4;
        el_type = 4;
        cells_nodes.reserve(tetra4_nodes.size());
        cells_nodes = tetra4_nodes;
        
        gauss_points_tetra4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
        n_gauss_cells = gauss_weight_cells.Length();
    }
    if (n_hexa8>0){
        n_cells = n_hexa8;
        el_type = 8;
        cells_nodes.reserve(hexa8_nodes.size());
        cells_nodes = hexa8_nodes;
        
        gauss_points_hexa8(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
        n_gauss_cells = gauss_weight_cells.Length();
    }
    if (n_tetra10>0){
        n_cells = n_tetra10;
        el_type= 10;
        cells_nodes.reserve(tetra10_nodes.size());
        cells_nodes = tetra10_nodes;
                
        gauss_points_tetra11(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
        n_gauss_cells = gauss_weight_cells.Length();
    }
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

void mesh::print_info(){
    std::cout << "n_cells" << "\t\t" << "el_type" << "\t\t" << "n_nodes" << "\t\t" << "n_faces" << "\n";
    std::cout << n_cells << "\t\t" << el_type << "\t\t" << n_nodes << "\t\t" << n_faces << "\n";
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
            local_dof_without_ghosts.push_back(3*j);
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
                    local_dof.push_back(3*node);
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
                    break;
                }
            case 4:
                if (mynpart[nodes[0]]==MyPID && mynpart[nodes[1]]==MyPID && mynpart[nodes[2]]==MyPID && mynpart[nodes[3]]==MyPID){
                    local_faces.push_back(i);
                    break;
                }
                break;
            case 6:
                if (mynpart[nodes[0]]==MyPID && mynpart[nodes[1]]==MyPID && mynpart[nodes[2]]==MyPID && mynpart[nodes[3]]==MyPID && mynpart[nodes[4]]==MyPID && mynpart[nodes[5]]==MyPID){
                    local_faces.push_back(i);
                    break;
                }
        }
    }
    n_local_faces = local_faces.size();
    }
}

void mesh::store_feinterp_tri(){
    
    int node, eglob;
    Epetra_SerialDenseVector N(face_type); //xi(3), eta(3)
    Epetra_SerialDenseMatrix JacobianMatrix(2,2), X(2,face_type), D(face_type,2);
    
    /*xi(0) = 1.0/6.0; xi(1) = 2.0/3.0; xi(2) = 1.0/6.0;
    eta(0) = 1.0/6.0; eta(1) = 1.0/6.0; eta(2) = 2.0/3.0;*/
    
    N_tri.Reshape(n_gauss_faces,face_type);
    D1_N_tri.Reshape(n_gauss_faces,face_type);
    D2_N_tri.Reshape(n_gauss_faces,face_type);
    detJac_tri.Reshape(n_local_faces,n_gauss_faces);
    
    switch (face_type){
        case 3:
            for (unsigned int gp=0; gp<n_gauss_faces; ++gp){
                tri3::shape_functions(N, xi_faces[gp], eta_faces[gp]);
                for (int inode=0; inode<face_type; ++inode){
                    N_tri(gp,inode) = N(inode);
                }
            }
            break;
        case 4:
            for (unsigned int gp=0; gp<n_gauss_faces; ++gp){
                quad4::shape_functions(N, xi_faces[gp], eta_faces[gp]);
                for (int inode=0; inode<face_type; ++inode){
                    N_tri(gp,inode) = N(inode);
                }
            }
            break;
        case 6:
            for (unsigned int gp=0; gp<n_gauss_faces; ++gp){
                tri6::shape_functions(N, xi_faces[gp], eta_faces[gp]);
                for (int inode=0; inode<face_type; ++inode){
                    N_tri(gp,inode) = N(inode);
                }
            }
            break;
    };
    
    for (unsigned int eloc=0; eloc<n_local_faces; ++eloc){
        eglob = local_faces[eloc];
        for (unsigned int inode=0; inode<face_type; inode++){
            node = faces_nodes[face_type*eglob+inode];
            X(0,inode) = nodes_coord[3*node];
            X(1,inode) = nodes_coord[3*node+1];
        }
        
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
            }
            for (unsigned int inode=0; inode<face_type; ++inode){
                D1_N_tri(gp,inode) = D(inode,0);
                D2_N_tri(gp,inode) = D(inode,1);
            }
            
            JacobianMatrix.Multiply('N','N',1.0,X,D,0.0);
            detJac_tri(eloc,gp) = fabs(JacobianMatrix(0,0)*JacobianMatrix(1,1) - JacobianMatrix(1,0)*JacobianMatrix(0,1));
        }
        
    }
}

void mesh::store_feinterp_tetra(){
    
    int node, eglob;
    double alpha, beta;
    Epetra_SerialDenseVector N(el_type); //, xi(4), eta(4), zeta(4);
    Epetra_SerialDenseMatrix JacobianMatrix(3,3), InverseJacobianMatrix(3,3), X(3,el_type), D(el_type,3), DX(el_type,3);
    
    /*switch (el_type){
        case 4:
            alpha = (5.0 - sqrt(5.0))/20.0;
            beta = (5.0 + 3.0*sqrt(5.0))/20.0;
            xi[0] = alpha; eta[0] = alpha; zeta[0] = alpha;
            xi[1] = alpha; eta[1] = alpha; zeta[1] = beta;
            xi[2] = alpha; eta[2] = beta; zeta[2] = alpha;
            xi[3] = beta; eta[3] = alpha; zeta[3] = alpha;
            break;
        case 8:
            alpha = std::sqrt(2.0/3.0);
            beta = std::sqrt(1.0/3.0);
            xi[0] = 0.0; eta[0] = alpha; zeta[0] = -beta;
            xi[1] = 0.0; eta[1] = -alpha; zeta[1] = -beta;
            xi[2] = alpha; eta[2] = 0.0; zeta[2] = beta;
            xi[3]= -alpha; eta[3] = 0.0; zeta[3] = beta;
            break;
        case 10:
            alpha = (5.0 - sqrt(5.0))/20.0;
            beta = (5.0 + 3.0*sqrt(5.0))/20.0;
            xi[0] = alpha; eta[0] = alpha; zeta[0] = alpha;
            xi[1] = alpha; eta[1] = alpha; zeta[1] = beta;
            xi[2] = alpha; eta[2] = beta; zeta[2] = alpha;
            xi[3] = beta; eta[3] = alpha; zeta[3] = alpha;
            break;
    }*/
    
    local_rows.Resize(3*el_type*n_local_cells);
    vol_tetra.Resize(n_local_cells);
    N_tetra.Reshape(el_type,n_gauss_cells);
    detJac_tetra.Reshape(n_local_cells,n_gauss_cells);
    DX_N_tetra.Reshape(n_gauss_cells*el_type,n_local_cells);
    DY_N_tetra.Reshape(n_gauss_cells*el_type,n_local_cells);
    DZ_N_tetra.Reshape(n_gauss_cells*el_type,n_local_cells);
    
    switch (el_type){
        case 4:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                tetra4::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_tetra(inode,gp) = N(inode);
                }
            }
            break;
        case 8:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                hexa8::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_tetra(inode,gp) = N(inode);
                }
            }
            break;
        case 10:
            for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
                tetra10::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
                for (int inode=0; inode<el_type; ++inode){
                    N_tetra(inode,gp) = N(inode);
                }
            }
            break;
    };
    
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
        
        vol_tetra(eloc) = 0.0;
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
            }
            jacobian_matrix(X,D,JacobianMatrix);
            std::cout << X << "\n";
            jacobian_det(JacobianMatrix,detJac_tetra(eloc,gp));
            dX_shape_functions(D,JacobianMatrix,detJac_tetra(eloc,gp),DX);
            vol_tetra(eloc) += gauss_weight_cells(gp)*detJac_tetra(eloc,gp);
            
            for (int inode=0; inode<el_type; ++inode){
                DX_N_tetra(gp+n_gauss_cells*inode,eloc) = DX(inode,0);
                DY_N_tetra(gp+n_gauss_cells*inode,eloc) = DX(inode,1);
                DZ_N_tetra(gp+n_gauss_cells*inode,eloc) = DX(inode,2);
            }
        }
    }
    
}


