#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "Epetra_SerialDenseSolver.h"
#include  <math.h>
#include "Compressible_Mooney_Transverse_Isotropic.hpp"
#include "Newton_Raphsonpp.hpp"

class objectiveFunction
{
private:
    
    Teuchos::ParameterList _paramList;
    Epetra_Comm * comm;
    Teuchos::RCP<Newton_Raphson> newton;
    Teuchos::RCP<NRL_ModelF> interface;
    
    unsigned int npoints;
    unsigned int nloads;
    std::vector<int> exp_cells;
    std::vector<double> my_xi;
    std::vector<double> my_eta;
    std::vector<double> my_exx;
    std::vector<double> my_eyy;
    std::vector<double> my_exy;
    std::vector<double> data_bc;
    
public:
        
    objectiveFunction(Epetra_Comm & Comm, Teuchos::ParameterList & paramList){
        
        comm = &Comm;
        interface = Teuchos::rcp(new NRL_ModelF(Comm,paramList));
        newton = Teuchos::rcp(new Newton_Raphson(*interface,paramList));
        
        std::string path_exp_points =
        Teuchos::getParameter<std::string>(paramList.sublist("Data"),"path_to_pts");
        std::string path_exp_deform =
        Teuchos::getParameter<std::string>(paramList.sublist("Data"),"path_to_def");
        
        if (Comm.MyPID()==0){
            std::cout << "npoints" << std::setw(30) << "nloads" << std::setw(30) << "local_npoints\n";
        }
        retrieve_data(path_exp_points, path_exp_deform);
        
        _paramList = paramList;
    }
    
    ~objectiveFunction(){
    }
    
    double value(Epetra_SerialDenseVector & x){
        
        double plyagl = 30.0*2.0*M_PI/360.0;
        interface->set_parameters(x);
        interface->set_plyagl(plyagl);
    
        double val = 0.0;
        newton->Initialization();
        for (unsigned int i=0; i<data_bc.size(); i=i+50){
            newton->setParameters(_paramList);
            newton->bc_disp=data_bc[i];
            int error = newton->Solve_with_Aztec(false);
            
            Epetra_SerialDenseVector exx_comp(exp_cells.size());
            Epetra_SerialDenseVector eyy_comp(exp_cells.size());
            Epetra_SerialDenseVector exy_comp(exp_cells.size());
            
            if (!error){
                compute_green_lagrange(*newton->x,exx_comp,eyy_comp,exy_comp);
            }
            else{
                if (comm->MyPID()==0){
                    std::cout << "Newton failed.\n";
                }
            }
            
            double absError = 0.0;
            double loadRef = 0.0;
            double partialLoadRef = 0.0;
            double partialAbsError = 0.0;
            for (unsigned int j=0; j<exp_cells.size(); ++j){
                partialAbsError += (exx_comp(j)*exx_comp(j)+eyy_comp(j)*eyy_comp(j)+2.0*exy_comp(j)*exy_comp(j) - my_exx[i+j*nloads]*my_exx[i+j*nloads]-my_eyy[i+j*nloads]*my_eyy[i+j*nloads]-2.0*my_exy[i+j*nloads]*my_exy[i+j*nloads]);
                partialLoadRef += my_exx[i+j*nloads]*my_exx[i+j*nloads]+my_eyy[i+j*nloads]*my_eyy[i+j*nloads]+2.0*my_exy[i+j*nloads]*my_exy[i+j*nloads];
            }
            comm->SumAll(&partialAbsError,&absError,1);
            comm->SumAll(&partialLoadRef,&loadRef,1);
            val += std::fabs(absError)/loadRef;
        }
        
        return val;
    }
    
    void retrieve_data(std::string & filepoints, std::string & filedefs){
        
        double testx, testy, testz, xi, eta, residual;
        unsigned int n_local_faces = interface->Mesh->n_local_faces;
        int node, result, e_gid;
        int nvert = interface->Mesh->face_type;
        
        Epetra_SerialDenseVector x(nvert), y(nvert), z(nvert);
        
        std::vector<double> data_xyz;
        std::vector<double> data_exx, data_eyy, data_exy;
        
        import_exp_points(filepoints, data_xyz);
        import_exp_def(filedefs, data_exx, data_eyy, data_exy);
        
        npoints = data_xyz.size()/3;
        nloads  = data_exx.size()/npoints;
        for (unsigned int p=0; p<npoints; ++p){
            testx = data_xyz[3*p+0]/1000.0;
            testy = data_xyz[3*p+1]/1000.0;
            testz = data_xyz[3*p+2]/1000.0;
            for (unsigned int e_lid=0; e_lid<n_local_faces; ++e_lid){
                e_gid = interface->Mesh->local_faces[e_lid];
                result = -1;
                for (unsigned int inode=0; inode<nvert; ++inode){
                    node = interface->Mesh->faces_nodes[nvert*e_gid+inode];
                    x(inode) = interface->Mesh->nodes_coord[3*node+0];
                    y(inode) = interface->Mesh->nodes_coord[3*node+1];
                    z(inode) = interface->Mesh->nodes_coord[3*node+2];
                }
                if (z(0)==testz && testx>=0.0 && testx<=50.0/1000.0 && testy>=0.0 && testy<=25.0/1000.0){
                    result = pnpoly(nvert,x,y,testx,testy);
                }
                if (result==1){
                    exp_cells.push_back(e_lid);
                    residual = inverse_isoparametric_mapping(testx,testy,x,y,xi,eta);
                    my_xi.push_back(xi);
                    my_eta.push_back(eta);
                    for (unsigned int load=0; load<nloads; ++load){
                        my_exx.push_back(data_exx[p+load*npoints]);
                        my_eyy.push_back(data_eyy[p+load*npoints]);
                        my_exy.push_back(data_exy[p+load*npoints]);
                    }
                }
            }
        }
        
        std::cout << npoints << std::setw(30) << nloads << std::setw(30) << exp_cells.size() << "\n";
    }
    
    void import_exp_points(std::string & filename, std::vector<double> & data_xyz){
        int ndata;
        double coordinate;
        std::ifstream file;
        
        file.open(filename);
        if (file.is_open()){
            file >> ndata;
            for (unsigned int i=0; i<3*ndata; ++i){
                file >> coordinate;
                data_xyz.push_back(coordinate);
            }
            file.close();
        }
        else{
            std::cout << "Couldn't open something\n";
        }
    }

    void import_exp_def(std::string & filename, std::vector<double> & data_exx, std::vector<double> & data_eyy, std::vector<double> & data_exy){
        
        int ndata, nload;
        double deformation, bc;
        std::ifstream file1, file2, file3, file4;
        
        file1.open(filename+"exx_id1.txt");
        file2.open(filename+"eyy_id1.txt");
        file3.open(filename+"exy_id1.txt");
        file4.open(filename+"bc_id1.txt");
        if (file1.is_open() && file2.is_open() && file3.is_open() && file4.is_open()){
            file1 >> ndata; file1 >> nload;
            file2 >> ndata; file2 >> nload;
            file3 >> ndata; file3 >> nload;
            for (unsigned int i=0; i<nload; ++i){
                for (unsigned int j=0; j<ndata; ++j){
                    file1 >> deformation;
                    data_exx.push_back(deformation);
                    file2 >> deformation;
                    data_eyy.push_back(deformation);
                    file3 >> deformation;
                    data_exy.push_back(deformation);
                }
                file4 >> bc;
                data_bc.push_back(bc);
            }
            file1.close();
            file2.close();
            file3.close();
            file4.close();
        }
        else{
            std::cout << "Couldn't open something with path " << filename << "\n";
        }
        
    }
    
    double inverse_isoparametric_mapping(double & testx, double & testy, Epetra_SerialDenseVector & x, Epetra_SerialDenseVector & y, double & xi, double & eta){
        
        Epetra_SerialDenseSolver solver;
        
        double rhs_inf = 1.0;
        int nvert = interface->Mesh->face_type;
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
    
    void compute_green_lagrange(Epetra_Vector & x, Epetra_SerialDenseVector & exx_comp, Epetra_SerialDenseVector & eyy_comp, Epetra_SerialDenseVector & exy_comp){
        
        Epetra_Vector u(*(interface->OverlapMap));
        u.Import(x, *(interface->ImportToOverlapMap), Insert);
        
        int e_gid, node;
        double det_jac;
        Epetra_SerialDenseMatrix matrix_X(2,interface->Mesh->face_type);
        Epetra_SerialDenseMatrix matrix_x(2,interface->Mesh->face_type);
        Epetra_SerialDenseMatrix D(interface->Mesh->face_type,2);
        Epetra_SerialDenseMatrix dx_shape_functions(interface->Mesh->face_type,2);
        Epetra_SerialDenseMatrix DX(interface->Mesh->face_type,2);
        Epetra_SerialDenseMatrix deformation_gradient(2,2);
        Epetra_SerialDenseMatrix JacobianMatrix(2,2);
        Epetra_SerialDenseMatrix InverseJacobianMatrix(2,2);
        Epetra_SerialDenseMatrix right_cauchy(2,2);
        
        for (unsigned e=0; e<exp_cells.size(); ++e){
            e_gid = interface->Mesh->local_faces[exp_cells[e]];
            
            for (unsigned int inode=0; inode<interface->Mesh->face_type; ++inode){
                node = interface->Mesh->faces_nodes[interface->Mesh->face_type*e_gid+inode];
                matrix_X(0,inode) = interface->Mesh->nodes_coord[3*node+0];
                matrix_X(1,inode) = interface->Mesh->nodes_coord[3*node+1];
                matrix_x(0,inode) = u[interface->OverlapMap->LID(3*node+0)] + interface->Mesh->nodes_coord[3*node+0];
                matrix_x(1,inode) = u[interface->OverlapMap->LID(3*node+1)] + interface->Mesh->nodes_coord[3*node+1];
            }
            
            switch (interface->Mesh->face_type){
                case 3:
                    tri3::d_shape_functions(D, my_xi[e], my_eta[e]);
                    break;
                case 4:
                    quad4::d_shape_functions(D, my_xi[e], my_eta[e]);
                    break;
                case 6:
                    tri6::d_shape_functions(D, my_xi[e], my_eta[e]);
                    break;
            }
            
            jacobian_matrix(matrix_X,D,JacobianMatrix);
            det_jac = fabs(JacobianMatrix(0,0)*JacobianMatrix(1,1) - JacobianMatrix(1,0)*JacobianMatrix(0,1));
            InverseJacobianMatrix(0,0) =  (1.0/det_jac)*JacobianMatrix(1,1);
            InverseJacobianMatrix(1,1) =  (1.0/det_jac)*JacobianMatrix(0,0);
            InverseJacobianMatrix(0,1) = -(1.0/det_jac)*JacobianMatrix(0,1);
            InverseJacobianMatrix(1,0) = -(1.0/det_jac)*JacobianMatrix(1,0);
            dx_shape_functions.Multiply('N','N',1.0,D,InverseJacobianMatrix,0.0);
            
            deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);
            right_cauchy.Multiply('T','N',1.0,deformation_gradient,deformation_gradient,0.0);
            
            exx_comp(e) = (1.0/2.0)*(right_cauchy(0,0)-1.0);
            eyy_comp(e) = (1.0/2.0)*(right_cauchy(1,1)-1.0);
            exy_comp(e) = (1.0/2.0)*right_cauchy(0,1);
        }
    }
    
};
