#include "distributenrldata.hpp"
#include <math.h>

distributenrldata::distributenrldata(mesh & Mesh){
    retrieve_data(Mesh);
}

distributenrldata::~distributenrldata(){
}

void distributenrldata::retrieve_data(mesh & Mesh){
    
    double testx, testy, testz, xi, eta, residual;
    unsigned int n_local_faces = Mesh.n_local_faces;
    int nvert = Mesh.face_type;
    int node, result, e_gid;
    
    Epetra_SerialDenseVector x(nvert), y(nvert), z(nvert);
    
    std::vector<double> data_xyz;
    std::vector<double> data_exx, data_eyy, data_exy;
    
    Teuchos::RCP<readnrldata> nrldata = Teuchos::rcp(new readnrldata());
    npoints = nrldata->npoints;
    nloads  = nrldata->nloads;
    boundaryconditions = nrldata->boundaryconditions;
    energy             = nrldata->energy;
    
    for (unsigned int p=0; p<npoints; ++p){
        testx = nrldata->points(p,0);
        testy = nrldata->points(p,1);
        testz = nrldata->points(p,2);
        for (unsigned int e_lid=0; e_lid<n_local_faces; ++e_lid){
            e_gid = Mesh.local_faces[e_lid];
            result = -1;
            for (unsigned int inode=0; inode<nvert; ++inode){
                node = Mesh.faces_nodes[nvert*e_gid+inode];
                x(inode) = Mesh.nodes_coord[3*node+0];
                y(inode) = Mesh.nodes_coord[3*node+1];
                z(inode) = Mesh.nodes_coord[3*node+2];
            }
            if (z(0)==testz && testx>=0.0 && testx<=50.0 && testy>=0.0 && testy<=25.0){
                result = pnpoly(nvert,x,y,testx,testy);
            }
            if (result==1){
                local_cells.push_back(e_lid);
                residual = inverse_isoparametric_mapping(testx,testy,x,y,xi,eta);
                local_xi.push_back(xi);
                local_eta.push_back(eta);
            }
        }
    }
}

double distributenrldata::inverse_isoparametric_mapping(double & testx, double & testy, Epetra_SerialDenseVector & x, Epetra_SerialDenseVector & y, double & xi, double & eta){
    
    Epetra_SerialDenseSolver solver;
    
    double rhs_inf = 1.0;
    int nvert = x.Length();;
    int iter = 0;
    
    Epetra_SerialDenseVector f(2);
    Epetra_SerialDenseVector dxi(2);
    Epetra_SerialDenseVector N(nvert);
    Epetra_SerialDenseMatrix D(nvert,2);
    Epetra_SerialDenseMatrix A(2,2);
    Epetra_SerialDenseMatrix mat_x(2,nvert);
    
    for (unsigned int i=0; i<nvert; ++i){
        mat_x(0,i) = x(i);
        mat_x(1,i) = y(i);
    }
    
    xi = 0.0; eta = 0.0;
    while (rhs_inf>1e-10){
        switch (nvert){
            case 3:
                tri3::shape_functions(N,xi,eta);
                tri3::d_shape_functions(D,xi,eta);
                break;
            case 4:
                quad4::shape_functions(N,xi,eta);
                quad4::d_shape_functions(D,xi,eta);
                break;
            case 6:
                tri6::shape_functions(N,xi,eta);
                tri6::d_shape_functions(D,xi,eta);
                break;
        }
        f(0) = testx;
        f(1) = testy;
        f.Multiply('N','N',-1.0,mat_x,N,1.0);
        
        iter++;
        if (iter>1){
            rhs_inf = f.NormInf();
            if (rhs_inf>1e8){
                std::cout << "Inverse Isoparametric Mapping: Newton-Raphson failed to converge.\n";
                break;
            }
        }
        if (iter>1000){
            std::cout << "Inverse Isoparametric Mapping: Iteration number exceeds 1000.\n";
            break;
        }
        
        A.Multiply('N','N',1.0,mat_x,D,0.0);
        
        solver.SetMatrix(A);
        solver.SetVectors(dxi,f);
        int error = solver.Solve();
        if (error){
            std::cout << "Inverse Isoparametric Mapping: Error with Epetra_SerialDenseSolver.\n";
        }
        xi += dxi(0);
        eta += dxi(1);
    }
    
    return rhs_inf;
}
