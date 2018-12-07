/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef NRL_PCA_LIKELIHOOD
#define NRL_PCA_LIKELIHOOD

#include  <math.h>
#include "compressible_Mooney_Transverse_Isotropic_Random_Field.hpp"
#include "newtonRaphson.hpp"
#include "distributenrldata.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>

class nrl_PCA_Likelihood
{
private:

    Teuchos::ParameterList      _paramList;
    Teuchos::RCP<newtonRaphson> newton;
    std::string                 fullOutputPath;

    Epetra_Comm              * comm;
    Epetra_Map               * MapExpPoints;
    Epetra_SerialDenseVector lb, ub;

public:

    Teuchos::RCP<tiMooneyRandomField> interface;
    Teuchos::RCP<distributenrldata>    nrldata;

    nrl_PCA_Likelihood(Epetra_Comm & Comm, Teuchos::ParameterList & paramList){
        comm                   = &Comm;
        _paramList             = paramList;
        std::string pathnrl    = Teuchos::getParameter<std::string>(paramList.sublist("nrldata"),"pathnrl");
        std::string station    = Teuchos::getParameter<std::string>(paramList.sublist("nrldata"),"station");
        std::string outputpath = Teuchos::getParameter<std::string>(paramList.sublist("nrldata"),"outputpath");
        fullOutputPath         = outputpath + station + "/";
        interface              = Teuchos::rcp(new tiMooneyRandomField(Comm,paramList));
        newton                 = Teuchos::rcp(new newtonRaphson(*interface,paramList));
        nrldata                = Teuchos::rcp(new distributenrldata(*interface->Mesh,pathnrl));

        MapExpPoints           = new Epetra_Map(-1,nrldata->global_id_faces.size(),&nrldata->global_id_faces[0],0,Comm);
    }

    nrl_PCA_Likelihood(){
    }

    int rnd(unsigned int                  nmc                  ,
            Epetra_IntSerialDenseVector & seeds                ,
            Epetra_SerialDenseVector    & mean_parameters      ,
            Epetra_SerialDenseVector    & exponents            ,
            Epetra_SerialDenseVector    & correlation_lengths  ,
            Epetra_SerialDenseVector    & coeff_of_variation   ,
            double                      & plyagl_deg           ,
            bool                          printNewtonIterations,
            bool                          printDisplacements   ,
            bool                          printDeformations    ,
            bool                          printRandomVariableY ,
            bool                          printDeformationsOnFaces)
    {
        Epetra_SerialDenseVector omega(6);
        omega(0) = coeff_of_variation(0);
        omega(1) = coeff_of_variation(1);
        omega(2) = coeff_of_variation(2);
        omega(3) = coeff_of_variation(3);
        omega(4) = correlation_lengths(0);
        omega(5) = correlation_lengths(1);

        double plyagl = plyagl_deg*2.0*M_PI/360.0;
        interface->setParameters(mean_parameters,exponents,omega);
        interface->set_plyagl(plyagl);

        if (comm->MyPID()==0){
          std::cout << "Generation of the Gaussian fields...";
        }
        interface->RandomFieldGenerator(seeds);
        if (comm->MyPID()==0){
          std::cout << "done.\n";
        }
        Epetra_Vector RandomVariableX(*MapExpPoints);
        RandomVariableX.PutScalar(0.0);

        Epetra_Vector DeformationsOnFaces_xx(*MapExpPoints);
        Epetra_Vector DeformationsOnFaces_yy(*MapExpPoints);
        Epetra_Vector DeformationsOnFaces_xy(*MapExpPoints);
        DeformationsOnFaces_xx.PutScalar(0.0);
        DeformationsOnFaces_yy.PutScalar(0.0);
        DeformationsOnFaces_xy.PutScalar(0.0);

        int error;
        newton->Initialization();
        for (unsigned int i=0; i<nrldata->boundaryconditions.Length(); ++i){
            newton->setParameters(_paramList);
            if(i==0){
                newton->bc_disp=nrldata->boundaryconditions(i);
            }
            else{
                newton->bc_disp=nrldata->boundaryconditions(i)-nrldata->boundaryconditions(i-1);
            }

            error = newton->Solve_with_Aztec(printNewtonIterations);

            if (!error){

                if (printDisplacements){
                  std::string path = fullOutputPath + "u_nmc=" + std::to_string(nmc) + "_angle=" + std::to_string(int(plyagl_deg)) + "_load=" + std::to_string(i) + "_delta=" + std::to_string(omega(0)) + "_l1" + std::to_string(omega(4)) + "_l2" + std::to_string(omega(5)) + ".mtx";
                  int flag = newton->print_newton_solution(path);
                  if (flag){
                    if (comm->MyPID()==0){
                      std::cout << "Failed printing displacements.";
                    }
                  }
                }

                if (printDeformations){
                  std::string path = fullOutputPath + "e_nmc" + std::to_string(nmc) + "_angle=" + std::to_string(int(plyagl_deg)) + "_load=" + std::to_string(i) + "_delta=" + std::to_string(omega(0)) + "_l1" + std::to_string(omega(4)) + "_l2" + std::to_string(omega(5)) + ".mtx";
                  double xi = 0.0;
                  int flag = interface->compute_green_lagrange(*newton->x,xi,xi,xi,path);
                  if (flag){
                    if (comm->MyPID()==0){
                      std::cout << "Failed printing deformations.";
                    }
                  }
                }

                Epetra_SerialDenseMatrix eij(nrldata->local_id_faces.size(),3);
                compute_green_lagrange(*newton->x,eij);

                if (printDeformationsOnFaces){
                  for (unsigned int i=0; i<nrldata->global_id_faces.size(); ++i){
                    DeformationsOnFaces_xx[MapExpPoints->LID(nrldata->global_id_faces[i])] = eij(i,0);
                    DeformationsOnFaces_yy[MapExpPoints->LID(nrldata->global_id_faces[i])] = eij(i,1);
                    DeformationsOnFaces_xy[MapExpPoints->LID(nrldata->global_id_faces[i])] = eij(i,2);
                  }
                  std::string pathDeformationsOnFaces_xx = fullOutputPath + "exx_nmc=" + std::to_string(nmc) + "_plyagl=" + std::to_string(int(plyagl_deg)) + "_s=" + std::to_string(i) + ".mtx";
                  std::string pathDeformationsOnFaces_yy = fullOutputPath + "eyy_nmc=" + std::to_string(nmc) + "_plyagl=" + std::to_string(int(plyagl_deg)) + "_s=" + std::to_string(i) + ".mtx";
                  std::string pathDeformationsOnFaces_xy = fullOutputPath + "exy_nmc=" + std::to_string(nmc) + "_plyagl=" + std::to_string(int(plyagl_deg)) + "_s=" + std::to_string(i) + ".mtx";
                  comm->Barrier();
                  int error1 = printIndicatorY(pathDeformationsOnFaces_xx, DeformationsOnFaces_xx);
                  int error2 = printIndicatorY(pathDeformationsOnFaces_yy, DeformationsOnFaces_yy);
                  int error3 = printIndicatorY(pathDeformationsOnFaces_xy, DeformationsOnFaces_xy);
                }

                for (unsigned int j=0; j<nrldata->local_id_faces.size(); ++j){
                    RandomVariableX[j] += eij(j,0)*eij(j,0)+eij(j,1)*eij(j,1)+2.0*eij(j,2)*eij(j,2);
                }
            }
            else{
                if (comm->MyPID()==0){
                    std::cout << "Newton failed at loading " << i << ".\n";
                    std::cout << "GIndicator(" << i << ") set to 0.0.\n";
                }
                error = 1;
                break;
            }
        }
        if (printRandomVariableY && !error){
          std::string pathToRandomVariableX = fullOutputPath + "RandomVariableY_angle=" + std::to_string(int(plyagl_deg)) + "_nmc=" + std::to_string(nmc) + ".mtx";
          int error = printIndicatorY(pathToRandomVariableX, RandomVariableX);
          return error;
        }
        else{
          return error;
        }
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

        for (unsigned int e=0; e<nrldata->local_id_faces.size(); ++e){
            e_gid = interface->Mesh->local_faces[nrldata->local_id_faces[e]];
            //if (interface->Mesh->faces_phases[e_gid]==92 || interface->Mesh->faces_phases[e_gid]==93){
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
        //}
    }

    /*void compute_force_reaction(Epetra_Vector & x){
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

      Epetra_SerialDenseVector TResult(2);
      Epetra_SerialDenseMatrix PiolaS(2,2);
      Epetra_SerialDenseVector TStress(2);
      Epetra_SerianDenseVector Normal(2);

      for (unsigned int e=0; e<interface->Mesh->local_faces.size(); ++e){
          e_gid = interface->Mesh->local_faces[e];
          if (interface->Mesh->faces_phases[e_gid]==90){
            for (unsigned int inode=0; inode<interface->Mesh->face_type; ++inode){
                node = interface->Mesh->faces_nodes[interface->Mesh->face_type*e_gid+inode];
                matrix_X(0,inode) = interface->Mesh->nodes_coord[3*node+0];
                matrix_X(1,inode) = interface->Mesh->nodes_coord[3*node+2];
                matrix_x(0,inode) = u[interface->OverlapMap->LID(3*node+0)] + interface->Mesh->nodes_coord[3*node+0];
                matrix_x(1,inode) = u[interface->OverlapMap->LID(3*node+2)] + interface->Mesh->nodes_coord[3*node+2];
            }
            for (unsigned int gp=0; gp<interface->Mesh->n_gauss_faces; ++gp){
              jacobian_matrix(matrix_X,D,JacobianMatrix);
              det_jac = fabs(JacobianMatrix(0,0)*JacobianMatrix(1,1) - JacobianMatrix(1,0)*JacobianMatrix(0,1));
              InverseJacobianMatrix(0,0) =  (1.0/det_jac)*JacobianMatrix(1,1);
              InverseJacobianMatrix(1,1) =  (1.0/det_jac)*JacobianMatrix(0,0);
              InverseJacobianMatrix(0,1) = -(1.0/det_jac)*JacobianMatrix(0,1);
              InverseJacobianMatrix(1,0) = -(1.0/det_jac)*JacobianMatrix(1,0);
              dx_shape_functions.Multiply('N','N',1.0,D,InverseJacobianMatrix,0.0);
              deformation_gradient.Multiply('N','N',1.0,matrix_x,dx_shape_functions,0.0);

              TStress.Multiplt('N','N',1.0,PiolaS,Normal,0.0);
              TResult.Multiply('N','N',interface->Mesh->gauss_weight_faces[gp],deformation_gradient,TStress,1.0);
            }
          }
      }
    }
    */

};
#endif
