#include "shinozukapp_2d.hpp"
#include <math.h>

shinozuka_2d::shinozuka_2d(int & nu, double & length_x, double & length_y):
psi_(0.0,1.0), phi_(0.0,1.0), order(nu), l1(length_x), l2(length_y){
}

shinozuka_2d::shinozuka_2d(){
}

template<typename typearg>
double shinozuka_2d::tau_beta(typearg & beta){
    double tau = -1.0 + (2.0*double(beta)-1.0)/double(order);
    return tau;
}

double shinozuka_2d::s_tau(double & tau){
    //double s = 0.0;
    //if (tau>=-1.0 && tau<=1.0){
    double s = (2.0/double(order))*(1.0-fabs(tau));
    //}
    return s;
}

void shinozuka_2d::generator(Epetra_Vector & v, mesh & Mesh){

    double c = std::cos(rotation);
    double s = std::sin(rotation);

    int node;
    double x, y;
    double ti, tj;
    double si, sj;
    double psi, phi, w, arg;
    v.PutScalar(0.0);

    for (int i=1; i<=order; ++i){
        ti = tau_beta<int>(i);
        si = s_tau(ti);
        for (int j=1; j<=order; ++j){
            tj = tau_beta<int>(j);
            sj = s_tau(tj);

                psi = psi_(rng);
                phi = phi_(rng);
                w   = std::sqrt(-std::log(psi));
                for (int inode=0; inode<v.MyLength(); ++inode){
                    node = Mesh.local_nodes_without_ghosts[inode];
                    x    = Mesh.nodes_coord[3*node+0];
                    y    = Mesh.nodes_coord[3*node+1];
                    arg  = 2.0*M_PI*phi + (M_PI/l1)*ti*x + (M_PI/l2)*tj*y;
                    //arg = 2.0*M_PI*phi + (M_PI/l1)*ti*(x*s+y*c) + (M_PI/l2)*tj*(x*c-y*s);
                    //arg = 2.0*M_PI*phi + (M_PI*ti*c/l1 + M_PI*tj*s/l2)*x + (-M_PI*ti*s/l1 + M_PI*tj*c/l2)*y;
                    v[inode] += std::sqrt(2.0*si*sj)*w*std::cos(arg);
                }

        }
    }

}

void shinozuka_2d::icdf_gamma(Epetra_Vector & V, Epetra_Vector & G, double & alpha, double & beta){
    for (unsigned int i=0; i<V.MyLength(); ++i){
        double erfx = boost::math::erf<double>(V[i]);
        double y = (1.0/2.0)*(1.0 + erfx);
        double yinv = boost::math::gamma_p_inv<double,double>(alpha,y);
        G[i] = yinv*beta;
    }
}

void shinozuka_2d::icdf_beta(Epetra_Vector & V, Epetra_Vector & B, double & tau1, double & tau2){
    for (unsigned int i=0; i<V.MyLength(); ++i){
        double erfx = boost::math::erf<double>(V[i]);
        double y = (1.0/2.0)*(1.0 + erfx);
        B[i] = boost::math::ibeta_inv<double,double,double>(tau1,tau2,y);
    }
}

shinozuka_2d::~shinozuka_2d(){
}
