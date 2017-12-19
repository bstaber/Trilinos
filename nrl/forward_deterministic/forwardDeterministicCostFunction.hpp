#ifndef FORWARDDETERMINISTICCOSTFUNCTION_HPP
#define FORWARDDETERMINISTICCOSTFUNCTION_HPP

#include  <math.h>
#include "Compressible_Mooney_Transverse_Isotropic.hpp"
#include "Newton_Raphsonpp.hpp"
#include "distributenrldata.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>

class forwardDeterministicCostFunction
{
private:
    
    Teuchos::ParameterList          _paramList;
    Epetra_Comm *                   comm;
    Teuchos::RCP<TIMooney>          interface;
    Teuchos::RCP<distributenrldata> nrldata;
    
public:
    
    forwardDeterministicCostFunction(Epetra_Comm & Comm, Teuchos::ParameterList & paramList){
        comm = &Comm;
        _paramList = paramList;
        std::string pathnrl = Teuchos::getParameter<std::string>(paramList.sublist("nrldata"),"pathnrl");
        interface  = Teuchos::rcp(new TIMooney(Comm,paramList));
        nrldata    = Teuchos::rcp(new distributenrldata(*interface->Mesh,pathnrl));
    }
    
    ~forwardDeterministicCostFunction(){
    }
    
    Epetra_SerialDenseVector value_id(Epetra_SerialDenseVector & x, int & id){
        double plyagl = nrldata->angles(id)*2.0*M_PI/360.0;
        interface->set_parameters(x);
        interface->set_plyagl(plyagl);
        
        Teuchos::RCP<Newton_Raphson> newton = Teuchos::rcp(new Newton_Raphson(*interface,_paramList));
        
        Epetra_SerialDenseVector val(nrldata->boundaryconditions.M());
        newton->Initialization();
        for (unsigned int i=0; i<nrldata->boundaryconditions.M(); ++i){
            newton->setParameters(_paramList);
            if(i==0){
                newton->bc_disp=nrldata->boundaryconditions(i,id);
            }
            else{
                newton->bc_disp=nrldata->boundaryconditions(i,id)-nrldata->boundaryconditions(i-1,id);
            }
            int error = newton->Solve_with_Aztec(false);
            
            Epetra_SerialDenseMatrix eij(nrldata->local_cells.size(),3);
            
            if (!error){
                get_green_lagrange(*newton->x,eij);
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
            val(i) = totalEnergy;
        }
        double xi = 0.0;
        std::string filename = "/home/s/staber/Trilinos_results/nrl/forward_deterministic/e22_id" + std::to_string(id) + ".mtx";
        interface->compute_green_lagrange(*newton->x, xi, xi, xi, filename);
        return val;
    }
    
    void get_green_lagrange(Epetra_Vector & x, Epetra_SerialDenseMatrix & eij){
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
