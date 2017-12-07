#include  <math.h>
#include "Compressible_Mooney_Transverse_Isotropic.hpp"
#include "Newton_Raphsonpp.hpp"
#include "distributenrldata.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>

class RandomSearch_DeterministicModel
{
private:
    
    Teuchos::ParameterList _paramList;
    Epetra_Comm * comm;
    Teuchos::RCP<Newton_Raphson> newton;
    Teuchos::RCP<TIMooney> interface;
    Teuchos::RCP<distributenrldata> nrldata;
    Epetra_SerialDenseVector lb, ub;
    
public:
    
    Epetra_SerialDenseVector solution;
        
    RandomSearch_DeterministicModel(Epetra_Comm & Comm, Teuchos::ParameterList & paramList){
        comm = &Comm;
        _paramList = paramList;
        interface  = Teuchos::rcp(new TIMooney(Comm,paramList));
        newton     = Teuchos::rcp(new Newton_Raphson(*interface,paramList));
        nrldata    = Teuchos::rcp(new distributenrldata(*interface->Mesh));
        
        solution.Resize(7); lb.Resize(7); ub.Resize(7);
        lb(0) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu1_inf");
        ub(0) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu1_sup");
        lb(1) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu2_inf");
        ub(1) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu2_sup");
        lb(2) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu3_inf");
        ub(2) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu3_sup");
        lb(3) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu4_inf");
        ub(3) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu4_sup");
        lb(4) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu5_inf");
        ub(4) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"mu5_sup");
        lb(5) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"beta4_inf");
        ub(5) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"beta4_sup");
        lb(6) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"beta5_inf");
        ub(6) = Teuchos::getParameter<double>(paramList.sublist("TIMooney"),"beta5_sup");
    }
    
    ~RandomSearch_DeterministicModel(){
    }
    
    double randomsearch(Epetra_SerialDenseVector & x, int & id, int & niter, double & tol){
        int eval = 1;
        double fval = value(x,id);
        printHeader();
        printStatus(eval,fval,x);
        
        int n = x.Length();
        Epetra_SerialDenseMatrix L(n,n);
        Epetra_SerialDenseVector v(n);
        
        boost::random::mt19937 rng(std::time(0));
        boost::random::normal_distribution<double> randn(0.0,1.0);
        boost::random::uniform_real_distribution<double> rand(0.0,1.0);
        
        double afval = fval;
        while(eval<=niter){
            for (unsigned int i=0; i<n; ++i){
                L(i,i) = 1.0;
            }
            v = x;
            while(fval>=afval){
                if (comm->MyPID()==0){
                    x = mvrandn(v,L,randn,rng);
                    for (unsigned int j=0; j<n; ++j){
                        if (x(j)<lb(j)){
                            x(j)=lb(j);
                        }
                        if (x(j)>ub(j)){
                            x(j)=ub(j);
                        }
                    }
                }
                comm->Broadcast(x.Values(),n,0);
                comm->Barrier();
                fval = value(x,id);
                printStatus(eval,fval,x);
                eval++;
            }
            afval = fval;
            if (fval<=tol){
                break;
            }
        }
        solution = x;
        return fval;
    }
    
    void printStatus(int eval, double value, Epetra_SerialDenseVector & x){
        comm->Barrier();
        if (comm->MyPID()==0){
            std::cout << eval << std::setw(15) << std::scientific << value << std::setw(15);
            for (unsigned int j=0; j<x.Length(); ++j){
                std::cout << std::setw(15) << x(j);
            }
            std::cout << "\n";
        }
    }
    
    void printHeader(){
        comm->Barrier();
        if (comm->MyPID()==0){
            std::cout << "Direct Random Search Algorithm\n";
            std::cout << "#eval" << std::setw(15) << "value";
            for (unsigned int i=0; i<8; ++i){
                std::cout << std::setw(15) << "x(" << i << ")";
            }
            std::cout << "\n";
        }
    }
    
    Epetra_SerialDenseVector mvrandn(Epetra_SerialDenseVector & x,
                                     Epetra_SerialDenseMatrix & L,
                                     boost::random::normal_distribution<double> & w,
                                     boost::random::mt19937 & rng){
        int n = x.Length();
        Epetra_SerialDenseVector z(n),g(n);
        for (unsigned int j=0; j<n; ++j){
            z(j) = w(rng);
        }
        g.Multiply('N','N',1.0,L,z,0.0);
        g+=x;
        return g;
    }
    
    double value(Epetra_SerialDenseVector & x, int & id){
        double plyagl = nrldata->angles(id)*2.0*M_PI/360.0;
        interface->set_parameters(x);
        interface->set_plyagl(plyagl);
        
        double val = 0.0;
        newton->Initialization();
        for (unsigned int i=0; i<nrldata->boundaryconditions.Length(); ++i){
            newton->setParameters(_paramList);
            if(i==0){
                newton->bc_disp=nrldata->boundaryconditions(i);
            }
            else{
                newton->bc_disp=nrldata->boundaryconditions(i)-nrldata->boundaryconditions(i-1);
            }
            int error = newton->Solve_with_Aztec(false);
            
            Epetra_SerialDenseMatrix eij(nrldata->local_cells.size(),3);
            
            if (!error){
                compute_green_lagrange(*newton->x,eij);
            }
            else{
                if (comm->MyPID()==0){
                    std::cout << "Newton failed.\n";
                }
            }
        
            double totalEnergy   = 0.0;
            double partialEnergy = 0.0;
            for (unsigned int j=0; j<nrldata->local_cells.size(); ++j){
                partialEnergy += eij(j,0)*eij(j,0)+eij(j,1)*eij(j,1)+2.0*eij(j,2)*eij(j,2);
            }
            comm->SumAll(&partialEnergy,&totalEnergy,1);
            val += totalEnergy;
        }
        val = val - nrldata->energy(id);
        val = fabs(val)/fabs(nrldata->energy(id));
        return val;
    }
    
    void compute_green_lagrange(Epetra_Vector & x, Epetra_SerialDenseMatrix & eij){
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
        
        for (unsigned e=0; e<nrldata->local_cells.size(); ++e){
            e_gid = interface->Mesh->local_faces[nrldata->local_cells[e]];
            for (unsigned int inode=0; inode<interface->Mesh->face_type; ++inode){
                node = interface->Mesh->faces_nodes[interface->Mesh->face_type*e_gid+inode];
                matrix_X(0,inode) = interface->Mesh->nodes_coord[3*node+0];
                matrix_X(1,inode) = interface->Mesh->nodes_coord[3*node+1];
                matrix_x(0,inode) = u[interface->OverlapMap->LID(3*node+0)] + interface->Mesh->nodes_coord[3*node+0];
                matrix_x(1,inode) = u[interface->OverlapMap->LID(3*node+1)] + interface->Mesh->nodes_coord[3*node+1];
            }
            switch (interface->Mesh->face_type){
                case 3:
                    tri3::d_shape_functions(D, nrldata->local_xi[e], nrldata->local_eta[e]);
                    break;
                case 4:
                    quad4::d_shape_functions(D, nrldata->local_xi[e], nrldata->local_eta[e]);
                    break;
                case 6:
                    tri6::d_shape_functions(D, nrldata->local_xi[e], nrldata->local_eta[e]);
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
            
            eij(e,0) = (1.0/2.0)*(right_cauchy(0,0)-1.0);
            eij(e,1) = (1.0/2.0)*(right_cauchy(1,1)-1.0);
            eij(e,2) = (1.0/2.0)*right_cauchy(0,1);
        }
    }
    
    
    
};
