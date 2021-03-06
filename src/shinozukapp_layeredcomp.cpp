/*
Brian Staber (brian.staber@gmail.com)
*/

#include "shinozukapp_layeredcomp.hpp"
#include <math.h>

shinozuka_layeredcomp::shinozuka_layeredcomp(int & nu, double & length_x, double & length_y, double & length_z):
psi_(0.0,1.0), phi_(0.0,1.0), order(nu), l1(length_x), l2(length_y), l3(length_z){
}

shinozuka_layeredcomp::shinozuka_layeredcomp(){
}

template<typename typearg>
double shinozuka_layeredcomp::tau_beta(typearg & beta){
    double tau = -1.0 + (2.0*double(beta)-1.0)/double(order);
    return tau;
}

double shinozuka_layeredcomp::s_tau(double & tau){
    double s = 0.0;
    if (tau>=-1.0 && tau<=1.0){
        s = (2.0/double(order))*(1.0-fabs(tau));
    }
    return s;
}

void shinozuka_layeredcomp::generator_gauss_points(Epetra_SerialDenseVector & v, mesh & Mesh, std::vector<int> & phase){

    int node, e_gid;
    int n_local_cells = Mesh.n_local_cells;
    int n_gauss_cells = Mesh.n_gauss_cells;
    double x, y, z;
    double xi, eta, zeta;
    double ti, tj, tk;
    double si, sj, sk;
    double psi, phi, w, arg;

    Epetra_SerialDenseVector vector_x(3);
    Epetra_SerialDenseVector shape_functions(Mesh.el_type);
    Epetra_SerialDenseMatrix matrix_X(3,Mesh.el_type);

    for (int il=0; il<32; ++il){
        for (int i=1; i<=order; ++i){
            ti = tau_beta<int>(i);
            si = s_tau(ti);
            for (int j=1; j<=order; ++j){
                tj = tau_beta<int>(j);
                sj = s_tau(tj);
                for (int k=1; k<=order; ++k){
                    tk = tau_beta<int>(k);
                    sk = s_tau(tk);
                    psi = psi_(rng);
                    phi = phi_(rng);
                    w = std::sqrt(-std::log(psi));
                    for (int e_lid=il; e_lid<n_local_cells; e_lid++){
                        e_gid = Mesh.local_cells[e_lid];
                        if(phase[e_gid]==il){
                            for (unsigned inode=0; inode<Mesh.el_type; ++inode){
                                node = Mesh.cells_nodes[Mesh.el_type*e_gid+inode];
                                matrix_X(0,inode) = Mesh.nodes_coord[3*node+0];
                                matrix_X(1,inode) = Mesh.nodes_coord[3*node+1];
                                matrix_X(2,inode) = Mesh.nodes_coord[3*node+2];
                            }
                            for (int gp=0; gp<n_gauss_cells; ++gp){
                                xi = Mesh.xi_cells(gp); eta = Mesh.eta_cells(gp); zeta = Mesh.zeta_cells(gp);
                                switch (Mesh.el_type){
                                    case 4:
                                        tetra4::shape_functions(shape_functions, xi, eta, zeta);
                                        break;
                                    case 8:
                                        hexa8::shape_functions(shape_functions, xi, eta, zeta);
                                        break;
                                    case 10:
                                        tetra10::shape_functions(shape_functions, xi, eta, zeta);
                                        break;
                                }
                                vector_x.Multiply('N','N',1.0,matrix_X,shape_functions,0.0);
                                arg = 2.0*M_PI*phi + (M_PI/l1)*ti*vector_x(0) + (M_PI/l2)*tj*vector_x(1) + (M_PI/l3)*tk*vector_x(2);
                                v(e_lid*n_gauss_cells+gp) += std::sqrt(2.0*si*sj*sk)*w*std::cos(arg);
                            }
                        }
                    }
                }
            }
        }
    }

}

void shinozuka_layeredcomp::generator_one_gauss_point(Epetra_SerialDenseVector & v, mesh & Mesh, std::vector<int> & phase, double & xi, double & eta, double & zeta){

    int node, e_gid;
    int n_local_cells = Mesh.n_local_cells;
    int n_gauss_cells = Mesh.n_gauss_cells;
    double x, y, z;
    double ti, tj, tk;
    double si, sj, sk;
    double psi, phi, w, arg;

    Epetra_SerialDenseVector vector_x(3);
    Epetra_SerialDenseVector shape_functions(Mesh.el_type);
    Epetra_SerialDenseMatrix matrix_X(3,Mesh.el_type);

    switch (Mesh.el_type){
        case 4:
            tetra4::shape_functions(shape_functions, xi, eta, zeta);
            break;
        case 8:
            hexa8::shape_functions(shape_functions, xi, eta, zeta);
            break;
        case 10:
            tetra10::shape_functions(shape_functions, xi, eta, zeta);
            break;
    }

    for (int il=0; il<32; ++il){
        for (int i=1; i<=order; ++i){
            ti = tau_beta<int>(i);
            si = s_tau(ti);
            for (int j=1; j<=order; ++j){
                tj = tau_beta<int>(j);
                sj = s_tau(tj);
                for (int k=1; k<=order; ++k){
                    tk = tau_beta<int>(k);
                    sk = s_tau(tk);
                    psi = psi_(rng);
                    phi = phi_(rng);
                    w = std::sqrt(-std::log(psi));
                    for (int e_lid=0; e_lid<n_local_cells; ++e_lid){
                        e_gid = Mesh.local_cells[e_lid];
                            if(phase[e_gid]==il){
                                for (unsigned inode=0; inode<Mesh.el_type; ++inode){
                                    node = Mesh.cells_nodes[Mesh.el_type*e_gid+inode];
                                    matrix_X(0,inode) = Mesh.nodes_coord[3*node+0];
                                    matrix_X(1,inode) = Mesh.nodes_coord[3*node+1];
                                    matrix_X(2,inode) = Mesh.nodes_coord[3*node+2];
                                }
                                vector_x.Multiply('N','N',1.0,matrix_X,shape_functions,0.0);
                                arg = 2.0*M_PI*phi + (M_PI/l1)*ti*vector_x(0) + (M_PI/l2)*tj*vector_x(1) + (M_PI/l3)*tk*vector_x(2);
                                v(e_lid) += std::sqrt(2.0*si*sj*sk)*w*std::cos(arg);
                            }
                    }
                }
            }
        }
    }

}

void shinozuka_layeredcomp::icdf_gamma(Epetra_Vector & V, Epetra_Vector & G, double & alpha, double & beta){
    for (unsigned int i=0; i<V.MyLength(); ++i){
        double erfx = boost::math::erf<double>(V[i]/std::sqrt(2.0));
        double y = (1.0/2.0)*(1.0 + erfx);
        double yinv = boost::math::gamma_p_inv<double,double>(alpha,y);
        G[i] = yinv*beta;
    }
}

void shinozuka_layeredcomp::icdf_beta(Epetra_Vector & V, Epetra_Vector & B, double & tau1, double & tau2){
    for (unsigned int i=0; i<V.MyLength(); ++i){
        double erfx = boost::math::erf<double>(V[i]/std::sqrt(2.0));
        double y = (1.0/2.0)*(1.0 + erfx);
        B[i] = boost::math::ibeta_inv<double,double,double>(tau1,tau2,y);
    }
}

shinozuka_layeredcomp::~shinozuka_layeredcomp(){
}
