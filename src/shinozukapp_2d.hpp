/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef SHINOZUKAPP_2D_HPP
#define SHINOZUKAPP_2D_HPP

#include "meshpp.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

class shinozuka_2d
{
public:

    shinozuka_2d();
    shinozuka_2d(int & nu, double & L1, double & L2);
    ~shinozuka_2d();

    template<typename typearg>
    double tau_beta(typearg & beta);
    double s_tau(double & tau);

    void generator(Epetra_Vector & v, mesh & Mesh);

    void icdf_gamma(Epetra_Vector & V, Epetra_Vector & G, double & alpha, double & beta);
    void icdf_beta(Epetra_Vector & V, Epetra_Vector & B, double & tau1, double & tau2);

    int order;
    double l1;
    double l2;
    double rotation = 0.0;
    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> phi_;
    boost::random::uniform_real_distribution<> psi_;
};
#endif
