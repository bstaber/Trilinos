#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "NRL_RandomField.hpp"
#include "Newton_Raphsonpp.hpp"
#include "Epetra_SerialDenseSolver.h"
#include  <math.h>

template<class Real, class XPrim=ROL::StdVector<Real>, class XDual=ROL::StdVector<Real> >
class objectiveFunction : public ROL::Objective<Real>{
    
//    typedef std::vector<Real> vector;
//    typedef Vector<Real>      V;
//    typedef typename vector::size_type uint;
    
private:
    
    template<class VectorType>
    Teuchos::RCP< const std::vector<Real> > getVector( const ROL::Vector<Real>& x ) {
        return Teuchos::dyn_cast<const VectorType>((x)).getVector();
    }

    template<class VectorType>
    Teuchos::RCP< std::vector<Real> > getVector( ROL::Vector<Real>& x ) {
        return Teuchos::dyn_cast<VectorType>(x).getVector();
    }
    
    Teuchos::RCP<Newton_Raphson> newton;
    
    std::vector<int> exp_cells;
    std::vector<double> my_xi;
    std::vector<double> my_eta;
    std::vector<double> my_exx;
    std::vector<double> my_eyy;
    std::vector<double> my_exy;
    
public:
    
    Teuchos::RCP<NRL_RandomFieldModel> interface;
    
    objectiveFunction(Epetra_Comm & Comm, Teuchos::ParameterList & paramList){
        
        interface = Teuchos::rcp(new NRL_RandomFieldModel(Comm,paramList));
        newton = Teuchos::rcp(new Newton_Raphson(*interface,paramList));
        
        std::string path_exp_points =
        Teuchos::getParameter<std::string>(paramList.sublist("Data"),"path_to_pts");
        std::string path_exp_deform =
        Teuchos::getParameter<std::string>(paramList.sublist("Data"),"path_to_def");
        
        retrieve_data(path_exp_points, path_exp_deform);
        
    }
    ~objectiveFunction(){
    }
    
    Real value(const ROL::Vector<Real> &x, Real &tol){
        
        Teuchos::RCP<const std::vector<Real> > xp = getVector<XPrim>(x);
        uint n = xp->size();
        
        double m1 = (*xp)[1];
        double m2 = (*xp)[2];
        double beta3 = (*xp)[3];
        double beta4 = (*xp)[4];
        double beta5 = (*xp)[5];
        double plyagl = 45.0*2.0*M_PI/360.0;

        newton->Initialization();
        int error = newton->Solve_with_Aztec();
        
        //compute deformation on the boundaries.
        Real val = 0.0;
        return val;
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
        double deformation;
        std::ifstream file1, file2, file3;
        
        file1.open(filename+"_exx.txt");
        file2.open(filename+"_eyy.txt");
        file3.open(filename+"_exy.txt");
        if (file1.is_open() && file2.is_open() && file3.is_open()){
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
            }
            file1.close();
            file2.close();
            file3.close();
        }
        else{
            std::cout << "Couldn't open something\n";
        }
        
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
            
        unsigned int npoints = data_xyz.size()/3;
        unsigned int nloads  = data_exx.size()/npoints;
        
        for (unsigned int p=0; p<npoints; ++p){
            testx = data_xyz[3*p+0];
            testy = data_xyz[3*p+1];
            testz = data_xyz[3*p+2];
            for (unsigned int e_lid=0; e_lid<n_local_faces; ++e_lid){
                e_gid = interface->Mesh->local_faces[e_lid];
                result = -1;
                for (unsigned int inode=0; inode<nvert; ++inode){
                    node = interface->Mesh->faces_nodes[nvert*e_gid+inode];
                    x(inode) = interface->Mesh->nodes_coord[3*node+0];
                    y(inode) = interface->Mesh->nodes_coord[3*node+1];
                    z(inode) = interface->Mesh->nodes_coord[3*node+2];
                }
                if (z(0)==testz){
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
    
    int pncell(Epetra_SerialDenseVector & p){
        
        int result, node, e_gid;
        int nvert = interface->Mesh->face_type;
        Epetra_SerialDenseVector x(nvert), y(nvert), z(nvert);
        double testx = p(0);
        double testy = p(1);
        double testz = p(2);
        
        for (unsigned int e_lid=0; e_lid<interface->Mesh->n_local_faces; ++e_lid){
            e_gid = interface->Mesh->local_faces[e_lid];
            for (unsigned int inode=0; inode<nvert; ++inode){
                node = interface->Mesh->faces_nodes[nvert*e_gid+inode];
                x(inode) = interface->Mesh->nodes_coord[3*node+0];
                y(inode) = interface->Mesh->nodes_coord[3*node+1];
                z(inode) = interface->Mesh->nodes_coord[3*node+2];
            }
            if (z(0)==testz){
                result = pnpoly(nvert,x,y,testx,testy);
            }
            if (result==1){
                e_gid = interface->Mesh->local_faces[e_lid];
                return e_gid;
            }
        }
        
        return -1;
    }
    
};
