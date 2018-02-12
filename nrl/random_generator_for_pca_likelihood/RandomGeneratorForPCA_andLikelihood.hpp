#ifndef RANDOMGENERATORFORPCA_ANDLIKELIHOOD
#define RANDOMGENERATORFORPCA_ANDLIKELIHOOD

#include  <math.h>
#include "Compressible_Mooney_Transverse_Isotropic_Random_Field.hpp"
#include "Newton_Raphsonpp.hpp"
#include "distributenrldata.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>

class RandomGeneratorForPCA_andLikelihood
{
private:

    Teuchos::ParameterList _paramList;
    Teuchos::RCP<Newton_Raphson> newton;

    Epetra_Comm * comm;
    Epetra_Map  * MapExpPoints;
    Epetra_SerialDenseVector lb, ub;

public:

    Teuchos::RCP<TIMooney_RandomField> interface;
    Teuchos::RCP<distributenrldata>    nrldata;

    RandomGeneratorForPCA_andLikelihood(Epetra_Comm & Comm, Teuchos::ParameterList & paramList){
        comm = &Comm;
        _paramList = paramList;
        std::string pathnrl = Teuchos::getParameter<std::string>(paramList.sublist("nrldata"),"pathnrl");
        interface  = Teuchos::rcp(new TIMooney_RandomField(Comm,paramList));
        newton     = Teuchos::rcp(new Newton_Raphson(*interface,paramList));
        nrldata    = Teuchos::rcp(new distributenrldata(*interface->Mesh,pathnrl));

        MapExpPoints = new Epetra_Map(-1,nrldata->local_cells.size(),&nrldata->local_cells[0],0,Comm);
    }

    RandomGeneratorForPCA_andLikelihood(){
    }

    Epetra_SerialDenseVector rnd(unsigned int                & nmc,
                                 Epetra_IntSerialDenseVector & seeds,
                                 Epetra_SerialDenseVector    & mean_parameters,
                                 Epetra_SerialDenseVector    & exponents,
                                 Epetra_SerialDenseVector    & correlation_lengths,
                                 Epetra_SerialDenseVector    & coeff_of_variation,
                                 double                      & plyagl,
                                 bool                          printNewtonIterations,
                                 bool                          printDisplacements,
                                 bool                          printDeformations)
    {
        Epetra_SerialDenseVector omega(6);
        omega(0) = coeff_of_variation(0);
        omega(1) = coeff_of_variation(1);
        omega(2) = coeff_of_variation(2);
        omega(3) = coeff_of_variation(3);
        omega(4) = correlation_lengths(0);
        omega(5) = correlation_lengths(1);

        interface->setParameters(mean_parameters,exponents,omega);
        interface->set_plyagl(plyagl);
        interface->RandomFieldGenerator(seeds);

        Epetra_Vector            GIndicatorY(*MapExpPoints);
        Epetra_SerialDenseVector GIndicatorZ(nrldata->boundaryconditions.Length());

        GIndicatorY.PutScalar(0.0);
        newton->Initialization();
        for (unsigned int i=0; i<nrldata->boundaryconditions.Length(); ++i){
            newton->setParameters(_paramList);
            if(i==0){
                newton->bc_disp=nrldata->boundaryconditions(i);
            }
            else{
                newton->bc_disp=nrldata->boundaryconditions(i)-nrldata->boundaryconditions(i-1);
            }

            int error = newton->Solve_with_Aztec(printNewtonIterations);

            if (!error){

                if (printDisplacements){
                  std::string path = "/home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likelihood/u_nmc=" + std::to_string(nmc) + "_angle=" + std::to_string(plyagl) + "_k=" + std::to_string(i) + ".mtx";
                  int flag = newton->print_newton_solution(path);
                  if (flag){
                    if (comm->MyPID()==0){
                      std::cout << "Failed printing displacements.";
                    }
                  }
                }

                if (printDeformations){
                  std::string path = "/home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likelihood/e_nmc" + std::to_string(nmc) + "_angle=" + std::to_string(plyagl) + "_k=" + std::to_string(i) + ".mtx";
                  double xi = 0.0;
                  int flag = interface->compute_green_lagrange(*newton->x,xi,xi,xi,path);
                  if (flag){
                    if (comm->MyPID()==0){
                      std::cout << "Failed printing deformations.";
                    }
                  }
                }

                Epetra_SerialDenseMatrix eij(nrldata->local_cells.size(),3);
                compute_green_lagrange(*newton->x,eij);

                double LIndicatorZ = 0.0;
                for (unsigned int j=0; j<nrldata->local_cells.size(); ++j){
                    GIndicatorY[j] += eij(j,0)*eij(j,0)+eij(j,1)*eij(j,1)+2.0*eij(j,2)*eij(j,2);
                    LIndicatorZ    += eij(j,0)*eij(j,0)+eij(j,1)*eij(j,1)+2.0*eij(j,2)*eij(j,2);
                }
                comm->SumAll(&LIndicatorZ,&GIndicatorZ(i),1);

            }
            else{
                if (comm->MyPID()==0){
                    std::cout << "Newton failed at loading " << i << ".\n";
                    std::cout << "GIndicator(" << i << ") set to 0.0.\n";
                }
            }

        }
        std::cout << GIndicatorY;
        int error = printIndicatorY("/home/s/staber/Trilinos_results/nrl/random_generator_for_pca_likelihood/IndicatorY_nmc=" + std::to_string(nmc) + ".mtx",GIndicatorY);
        return GIndicatorZ;
    }

    int printIndicatorY(std::string filename, Epetra_Vector & X){
      int NumTargetElements = 0;
      if (comm->MyPID()==0){
          NumTargetElements = nrldata->npoints;
      }
      Epetra_Map         MapOnRoot(-1,NumTargetElements,0,*comm);
      Epetra_Export      ExportOnRoot(*MapExpPoints,MapOnRoot);
      Epetra_MultiVector lhs_root(MapOnRoot,true);
      lhs_root.Export(X,ExportOnRoot,Insert);

      int error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
      return error;
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
        Epetra_SerialDenseMatrix deformation_gradient(2,2);
        Epetra_SerialDenseMatrix JacobianMatrix(2,2);
        Epetra_SerialDenseMatrix InverseJacobianMatrix(2,2);
        Epetra_SerialDenseMatrix right_cauchy(2,2);

        for (unsigned int e=0; e<nrldata->local_cells.size(); ++e){
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
