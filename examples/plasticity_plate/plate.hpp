/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef plate_HPP
#define plate_HPP

#include "plasticitySmallStrains.hpp"

class plate : public plasticitySmallStrains
{
public:

    double E = 210.0;
    double nu = 0.3;
    double lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
    double mu = E/(2.0*(1.0+nu));
    double kappa = lambda + 2.0*mu/3.0;

    double R0 = 1000.0;
    double H  = 10000.0;

    Epetra_SerialDenseMatrix K4;
    Epetra_SerialDenseMatrix J4;

    Teuchos::ParameterList * Krylov;

    plate(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){

      plasticitySmallStrains::initialize(comm, Parameters);

      K4.Reshape(6,6);
      J4.Reshape(6,6);

      double c = 1.0/3.0;
      double d = 2.0/3.0;

      J4(0,0) = c;   J4(0,1) = c;   J4(0,2) = c;   J4(0,3) = 0.0; J4(0,4) = 0.0; J4(0,5) = 0.0;
      J4(1,0) = c;   J4(1,1) = c;   J4(1,2) = c;   J4(1,3) = 0.0; J4(1,4) = 0.0; J4(1,5) = 0.0;
      J4(2,0) = c;   J4(2,1) = c;   J4(2,2) = c;   J4(2,3) = 0.0; J4(2,4) = 0.0; J4(2,5) = 0.0;
      J4(3,0) = 0.0; J4(3,1) = 0.0; J4(3,2) = 0.0; J4(3,3) = 0.0; J4(3,4) = 0.0; J4(3,5) = 0.0;
      J4(4,0) = 0.0; J4(4,1) = 0.0; J4(4,2) = 0.0; J4(4,3) = 0.0; J4(4,4) = 0.0; J4(4,5) = 0.0;
      J4(5,0) = 0.0; J4(5,1) = 0.0; J4(5,2) = 0.0; J4(5,3) = 0.0; J4(5,4) = 0.0; J4(5,5) = 0.0;

      K4(0,0) = d;    K4(0,1) = -c;   K4(0,2) = -c;  K4(0,3) = 0.0; K4(0,4) = 0.0; K4(0,5) = 0.0;
      K4(1,0) = -c;   K4(1,1) = d;    K4(1,2) = -c;  K4(1,3) = 0.0; K4(1,4) = 0.0; K4(1,5) = 0.0;
      K4(2,0) = -c;   K4(2,1) = -c;   K4(2,2) = d;   K4(2,3) = 0.0; K4(2,4) = 0.0; K4(2,5) = 0.0;
      K4(3,0) = 0.0;  K4(3,1) = 0.0;  K4(3,2) = 0.0; K4(3,3) = 1.0; K4(3,4) = 0.0; K4(3,5) = 0.0;
      K4(4,0) = 0.0;  K4(4,1) = 0.0;  K4(4,2) = 0.0; K4(4,3) = 0.0; K4(4,4) = 1.0; K4(4,5) = 0.0;
      K4(5,0) = 0.0;  K4(5,1) = 0.0;  K4(5,2) = 0.0; K4(5,3) = 0.0; K4(5,4) = 0.0; K4(5,5) = 1.0;

    }

    ~plate(){
    }

    void constitutive_problem(const unsigned int & elid, const unsigned int & igp,
                              const Epetra_SerialDenseVector & DETO, Epetra_SerialDenseVector & SIG,
                              double & EPCUM, Epetra_SerialDenseMatrix & TGM){

      Epetra_SerialDenseVector SIGTR(6);
      Epetra_SerialDenseVector DSIG(6);
      Epetra_SerialDenseVector DEVTR(6);
      Epetra_SerialDenseVector EYE(6);
      EYE(0) = 1.0; EYE(1) = 1.0; EYE(2) = 1.0; EYE(3) = 0.0; EYE(4) = 0.0; EYE(5) = 0.0;

      DSIG.Multiply('N','N',1.0,ELASTICITY,DETO,0.0);
      for (unsigned int k=0; k<6; ++k) SIGTR(k) = SIG(k) + DSIG(k);

      double pressure = (1.0/3.0)*(SIGTR(0)+SIGTR(1)+SIGTR(2));

      DEVTR = SIGTR;
      DEVTR(0) -= pressure;
      DEVTR(1) -= pressure;
      DEVTR(2) -= pressure;

      double devdev = DEVTR(0)*DEVTR(0) + DEVTR(1)*DEVTR(1) + DEVTR(2)*DEVTR(2)
                    + DEVTR(3)*DEVTR(3) + DEVTR(4)*DEVTR(4) + DEVTR(5)*DEVTR(5);
      double qtrial = std::sqrt(1.5*devdev);
      double yield = qtrial - R0 - H*EPCUM;

      Epetra_SerialDenseMatrix NN(6,6);
      for (unsigned int i=0; i<6; ++i) {
        for (unsigned int j=0; j<6; ++j) {
          NN(i,j) = DEVTR(i)*DEVTR(j)/devdev;
        }
      }

      SIG = SIGTR;
      get_elasticity_tensor(elid, igp, TGM);
      if (yield>1.0e-10) {
        double gamma = yield/(3.0*mu+H);
        double beta  = 3.0*mu*gamma/qtrial;
        double alpha = 1.0-beta;

        for (unsigned int k=0; k<6; ++k) SIG(k) -= beta*DEVTR(k);

        EPCUM += gamma;

        for (unsigned int i=0; i<6; ++i) {
          for (unsigned int j=0; j<6; ++j) {
            TGM(i,j) =  + 3.0*kappa+J4(i,j) + 2.0*mu*alpha*K4(i,j) + 6.0*mu*mu*(gamma/qtrial - (1.0/(3.0*mu+H)))*NN(i,j);
          }
        }
      }

    }

    void get_elasticity_tensor(const unsigned int & e_lid, const unsigned int & gp, Epetra_SerialDenseMatrix & tgm){
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
        int dof = 1;
        double coord;
        unsigned int node;
        for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
            node = Mesh->local_nodes[i];
            coord = Mesh->nodes_coord[3*node+dof];
            if(coord>=0.0-1.0e-6 & coord<=0.0+1.0e-6){
                n_bc_dof+=3;
            }
            if(coord>=10.0-1.0e-6 & coord<=10.0+1.0e-6){
                n_bc_dof+=1;
            }
        }

        int indbc = 0;
        dof_on_boundary = new int [n_bc_dof];
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+dof];
            if (coord>=0.0-1.0e-6 & coord<=0.0+1.0e-6){
                dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+1] = 3*inode+1;
                dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=3;
            }
            if (coord>=10.0-1.0e-6 & coord<=10.0+1.0e-6){
                //dof_on_boundary[indbc+0] = 3*inode+0;
                dof_on_boundary[indbc+0] = 3*inode+1;
                //dof_on_boundary[indbc+2] = 3*inode+2;
                indbc+=1;
            }
        }
    }

    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        Epetra_MultiVector v(*StandardMap,true);
        v.PutScalar(0.0);

        int node;
        int dof = 1;
        double coord;
        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+dof];
            if (coord==10.0){
                v[0][StandardMap->LID(3*node+dof)] = displacement;
            }
        }

        Epetra_MultiVector rhs_dir(*StandardMap,true);
        K.Apply(v,rhs_dir);
        F.Update(-1.0,rhs_dir,1.0);

        for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
            node = Mesh->local_nodes[inode];
            coord = Mesh->nodes_coord[3*node+dof];
            if (coord>=0.0-1.0e-6 & coord<=0.0+1.0e-6){
                F[0][StandardMap->LID(3*node+0)] = 0.0;
                F[0][StandardMap->LID(3*node+1)] = 0.0;
                F[0][StandardMap->LID(3*node+2)] = 0.0;
            }
            if (coord>=10.0-1.0e-6 & coord<=10.0+1.0e-6){
                //F[0][StandardMap->LID(3*node+0)]   = 0.0;
                F[0][StandardMap->LID(3*node+dof)] = displacement;
                //F[0][StandardMap->LID(3*node+2)]   = 0.0;
            }
        }
        ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
    }
};
#endif
