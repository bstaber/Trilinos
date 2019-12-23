/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef LINEARPATCHTEST_HPP
#define LINEARPATCHTEST_HPP

#include "plasticitySmallStrains.hpp"

class linearPatchTest : public plasticitySmallStrains
{
public:

    double E = 210000.0;
    double nu = 0.3;
    double lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
    double mu = E/(2.0*(1.0+nu));
    double kappa = lambda + 2.0*mu/3.0;

    double R0 = 1000.0;
    double H  = 10000.0;

    Teuchos::ParameterList * Krylov;

    linearPatchTest(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){

    plasticitySmallStrains::initialize(comm, Parameters);

    }

    ~linearPatchTest(){
    }

    void constitutive_problem(const unsigned int & elid, const unsigned int & igp,
                              Epetra_SerialDenseVector & EEL, Epetra_SerialDenseVector & SIG,
                              Epetra_SerialDenseMatrix & TGM){
      //int gid = Mesh->local_cells[elid];
      //int idx = int(gid*Mesh->n_gauss_cells+igp)
      //GaussMap.LID(id)
      // implement the elastic response just for the first test
      TGM = ELASTICITY;
      SIG.Multiply('N','N',1.0,ELASTICITY,EEL,0.0);

      Epetra_SerialDenseVector DEV(6);
      Epetra_SerialDenseVector EYE(6);
      EYE(0) = 1.0; EYE(1) = 1.0; EYE(2) = 1.0; EYE(3) = 0.0; EYE(4) = 0.0; EYE(5) = 0.0;
      double pressure = (1.0/3.0)*(SIG(0)+SIG(1)+SIG(2));
      DEV = SIG;
      DEV(0) -= pressure;
      DEV(1) -= pressure;
      DEV(2) -= pressure;

      double qtrial = std::sqrt((3.0/2.0)*(DEV(0)*DEV(0) + DEV(1)*DEV(1) + DEV(2)*DEV(2) + DEV(3)*DEV(3)
                                         + DEV(4)*DEV(4) + DEV(5)*DEV(5)));

      double yield = qtrial - R0;
      if (yield>1.0e-10) {
        //std::cout << "yield = " << yield << std::endl;
        double gamma = yield/(3.0*mu+H);
        double alpha = (1.0 - (3.0*mu*gamma/qtrial));
        for (unsigned int k=0; k<6; ++k) {
          DEV(k) *= alpha;
          SIG(k) = DEV(k) + pressure*EYE(k);
          EEL(k) = (1.0/(2.0*mu))*DEV(k) + (1.0/(3.0*kappa))*pressure*EYE(k);
        }
      }
    }

    void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tgm){
        int e_gid = Mesh->local_cells[e_lid];
        int n_gauss_cells = Mesh->n_gauss_cells;

        double c1 = lambda+2.0*mu;

        tgm(0,0) = c1;     tgm(0,1) = lambda; tgm(0,2) = lambda; tgm(0,3) = 0.0;     tgm(0,4) = 0.0;    tgm(0,5) = 0.0;
        tgm(1,0) = lambda; tgm(1,1) = c1;     tgm(1,2) = lambda; tgm(1,3) = 0.0;     tgm(1,4) = 0.0;    tgm(1,5) = 0.0;
        tgm(2,0) = lambda; tgm(2,1) = lambda; tgm(2,2) = c1;     tgm(2,3) = 0.0;     tgm(2,4) = 0.0;    tgm(2,5) = 0.0;
        tgm(3,0) = 0.0;    tgm(3,1) = 0.0;    tgm(3,2) = 0.0;    tgm(3,3) = 2.0*mu;  tgm(3,4) = 0.0;    tgm(3,5) = 0.0;
        tgm(4,0) = 0.0;    tgm(4,1) = 0.0;    tgm(4,2) = 0.0;    tgm(4,3) = 0.0;     tgm(4,4) = 2.0*mu; tgm(4,5) = 0.0;
        tgm(5,0) = 0.0;    tgm(5,1) = 0.0;    tgm(5,2) = 0.0;    tgm(5,3) = 0.0;     tgm(5,4) = 0.0;    tgm(5,5) = 2.0*mu;

    }

    void setup_dirichlet_conditions(){
        n_bc_dof = 0;
        double x,y,z;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            //if(x==0.0||y==0.0||z==0.0||x==1.0||y==1.0||z==1.0){
            if(x==0.0||x==1.0){
                n_bc_dof+=3;
            }
        }

        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            //if(x==0.0||y==0.0||z==0.0||x==1.0||y==1.0||z==1.0){
            if(x==0.0||x==1.0){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
        }
    }

    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & TIME){
        Epetra_MultiVector v(*StandardMap,true);
        v.PutScalar(0.0);

        int node;
        double x,y,z;
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            //if (x==0.0||y==0.0||z==0.0||x==1.0||y==1.0||z==1.0){
            if(x==0.0||x==1.0){
                Epetra_SerialDenseVector u(3);
                u = exactSolution(x,y,z);
                v[0][StandardMap->LID(3*node+0)] = TIME*u(0); //0.1*(0.1*x + 0.08*y + 0.05*z + 0.04*x*y + 0.03*y*z + 0.08*z*x);
                v[0][StandardMap->LID(3*node+1)] = TIME*u(1); //x*y; //y*x; //0.1*(0.05*x + 0.04*y + 0.1*z + 0.07*x*y + 0.03*y*z + 0.08*z*x);
                v[0][StandardMap->LID(3*node+2)] = TIME*u(2); //x*y*z; //x*y*z; //0.1*(0.75*x + 0.09*y + 0.06*z + 0.07*x*y + 0.03*y*z + 0.08*z*x);
            }
        }

        Epetra_MultiVector rhs_dir(*StandardMap,true);
        K.Apply(v,rhs_dir);
        F.Update(-1.0,rhs_dir,1.0);

        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            x = Mesh->nodes_coord[3*node+0];
            y = Mesh->nodes_coord[3*node+1];
            z = Mesh->nodes_coord[3*node+2];
            //if (x==0.0||y==0.0||z==0.0||x==1.0||y==1.0||z==1.0){
            if(x==0.0||x==1.0){
                F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
                F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
                F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }

    double errorL2(Epetra_Vector & uStandardMap){

        Epetra_Vector u(*OverlapMap);
        u.Import(uStandardMap, *ImportToOverlapMap, Insert);
        double totalError;
        double error = 0.0;
        double normVH;
        double gauss_weight;
        int n_gauss_points = Mesh->n_gauss_cells;
        int e_gid, node;

        //Epetra_SerialDenseVector epsilon(6);
        //Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
        Epetra_SerialDenseVector uExact(3), vH(3);
        Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3), X_I(3,Mesh->el_type), u_I(3,Mesh->el_type);
        Epetra_SerialDenseMatrix u_G(3,n_gauss_points), x_G(3,n_gauss_points);

        for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
            e_gid = Mesh->local_cells[e_lid];
            for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
                X_I(0,inode) = Mesh->nodes_coord[3*node+0];
                X_I(1,inode) = Mesh->nodes_coord[3*node+1];
                X_I(2,inode) = Mesh->nodes_coord[3*node+2];
                u_I(0,inode) = u[OverlapMap->LID(3*node+0)];
                u_I(1,inode) = u[OverlapMap->LID(3*node+1)];
                u_I(2,inode) = u[OverlapMap->LID(3*node+2)];
            }
            x_G.Multiply('N','N',1.0,X_I,Mesh->N_cells,0.0);
            u_G.Multiply('N','N',1.0,u_I,Mesh->N_cells,0.0);
            for (unsigned int gp=0; gp<n_gauss_points; ++gp){
                gauss_weight = Mesh->gauss_weight_cells(gp);
                uExact = exactSolution(x_G(0,gp),x_G(1,gp),x_G(2,gp));
                vH(0) = uExact(0) - u_G(0,gp);
                vH(1) = uExact(1) - u_G(1,gp);
                vH(2) = uExact(2) - u_G(2,gp);
                normVH = vH.Norm2();
                /*for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
                 dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,e_lid);
                 dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,e_lid);
                 dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,e_lid);
                 }
                 compute_B_matrices(dx_shape_functions,matrix_B);
                 epsilon.Multiply('N','N',1.0,matrix_B,vector_u,0.0);*/
                error += gauss_weight*normVH*normVH*Mesh->detJac_cells(e_lid,gp);
            }
        }
        Comm->SumAll(&error,&totalError,1);
        totalError = std::sqrt(totalError);
        return totalError;
    }

    Epetra_SerialDenseVector exactSolution(double & x1, double & x2, double & x3){
        Epetra_SerialDenseVector u(3);
        u(0) = x1; //0.1*x1+0.2*x2+0.4*x3;
        u(1) = 0.0; //0.4*x1+0.5*x2+0.1*x3;
        u(2) = 0.0; //0.05*x1+0.25*x2+0.65*x3;
        return u;
    }
};


#endif
