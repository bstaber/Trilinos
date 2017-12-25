#ifndef RANDOMSEARCH_DETERMINSITICMODEL_HPP
#define RANDOMSEARCH_DETERMINSITICMODEL_HPP

#include <math.h>
#include <random>
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

    Teuchos::ParameterList          _paramList;
    Epetra_Comm *                   comm;
    Teuchos::RCP<Newton_Raphson>    newton;
    Teuchos::RCP<TIMooney>          interface;
    Teuchos::RCP<distributenrldata> nrldata;

public:

    Epetra_SerialDenseVector solution;

    RandomSearch_DeterministicModel(Epetra_Comm & Comm, Teuchos::ParameterList & paramList){
        comm = &Comm;
        _paramList = paramList;
        std::string pathnrl = Teuchos::getParameter<std::string>(paramList.sublist("nrldata"),"pathnrl");
        interface  = Teuchos::rcp(new TIMooney(Comm,paramList));
        newton     = Teuchos::rcp(new Newton_Raphson(*interface,paramList));
        nrldata    = Teuchos::rcp(new distributenrldata(*interface->Mesh,pathnrl));

        solution.Resize(7);
    }

    ~RandomSearch_DeterministicModel(){
    }

    double randomsearch(Epetra_SerialDenseVector & x, int & niter, double & tol){
        int eval = 1;
        double fval = value(x);
        printHeader();
        printStatus(eval,fval,x);

        int n = x.Length();
        Epetra_SerialDenseMatrix L(n,n);
        Epetra_SerialDenseVector v(n), lb(n), ub(n);

        for (unsigned int i=0; i<n; ++i){
            lb(i) = 0.10*x(i);
            ub(i) = 1.90*x(i);
        }

        static std::random_device rd;
        boost::random::mt19937 rng(rd());
        boost::random::normal_distribution<double>       randn(0.0,1.0);
        boost::random::uniform_real_distribution<double> rand(0.0,1.0);

        double fval, test;
        for (unsigned int eval=0; eval<niter; ++eval){

            for (unsigned int i=0; i<n; ++i){
                L(i,i) = 0.10*x(i);
            }
            v = x;

            if (comm->MyPID()==0){
              int flag = 1;
              while (flag==1){
                flag = 0;
                x = mvrandn(v,L,randn,rng);
                for (unsigned int j=0; j<n; ++j){
                  if (x(j)<lb(j) || x(j)>ub(j)){
                    flag = 1;
                    break;
                  }
                }
              }
            }
            comm->Broadcast(x.Values(),n,0);
            comm->Barrier();
            test = value(x);
            if (test<fval){
              printStatus(eval,fval,x);
              fval      = test;
              solultion = x;
            }
          }
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
            std::cout << "Direct Random Search Algorithm: " << nrldata->boundaryconditions.M() << " loads and " << nrldata->angles.Length() << " experimental tests.\n";
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

    double value(Epetra_SerialDenseVector & x){
        double vali = 0.0;
        double val  = 0.0;
        Epetra_SerialDenseVector angles(4);
        Epetra_IntSerialDenseMatrix angleToID(4,2);
        angles(0) = 2.0*M_PI*15.0/360.0;
        angles(1) = 2.0*M_PI*30.0/360.0;
        angles(2) = 2.0*M_PI*60.0/360.0;
        angles(3) = 2.0*M_PI*75.0/360.0;
        angleToID(0,0) = 4; angleToID(0,1) = 5;
        angleToID(1,0) = 0; angleToID(1,1) = 3;
        angleToID(2,0) = 1; angleToID(2,1) = 2;
        angleToID(3,0) = 6; angleToID(3,1) = 7;

        for (int i=0; i<4; ++i){
            vali = value_angle(x,angles(i),angleToID(i,0),angleToID(i,1));
            val += vali;
        }
        return val;
    }

    double value_angle(Epetra_SerialDenseVector & x, double & plyagl, int & id1, int & id2){
        //double plyagl = nrldata->angles(id)*2.0*M_PI/360.0;
        interface->set_parameters(x);
        interface->set_plyagl(plyagl);

        Epetra_SerialDenseVector meanEnergy(nrldata->boundaryconditions.Length());
        for (unsigned int i=0; i<meanEnergy.Length(); ++i){
            meanEnergy(i) = 0.5*(nrldata->energy(i,id1) + nrldata->energy(i,id2));
        }

        double val    = 0.0;
        double valref = 0.0;
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
            val    += (totalEnergy-meanEnergy(i))*(totalEnergy-meanEnergy(i));
            valref += meanEnergy(i)*meanEnergy(i);
        }
        val = val/valref;
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
#endif
