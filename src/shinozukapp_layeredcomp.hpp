#ifndef SHINOZUKAPP_LAYEREDCOMP_HPP
#define SHINOZUKAPP_LAYEREDCOMP_HPP

#include "meshpp.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

class shinozuka_layeredcomp
{
public:
    
    shinozuka_layeredcomp();
    shinozuka_layeredcomp(int & nu, double & L1, double & L2, double & L3);
    ~shinozuka_layeredcomp();
    
    template<typename typearg>
    double tau_beta(typearg & beta);
    double s_tau(double & tau);
    
    void generator_gauss_points(Epetra_SerialDenseVector & v, mesh & Mesh);
    void generator_one_gauss_point(Epetra_SerialDenseVector & v, mesh & Mesh, double & xi, double & eta, double & zeta);
    
    void icdf_gamma(Epetra_Vector & V, Epetra_Vector & G, double & alpha, double & beta);
    void icdf_beta(Epetra_Vector & V, Epetra_Vector & B, double & tau1, double & tau2);
        
    int order;
    double l1;
    double l2;
    double l3;
    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> phi_;
    boost::random::uniform_real_distribution<> psi_;
};
#endif
