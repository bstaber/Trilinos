#ifndef SHINOZUKAPP_LAYEREDCOMP_2D_HPP
#define SHINOZUKAPP_LAYEREDCOMP_2D_HPP

#include "meshpp.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

class shinozuka_layeredcomp_2d
{
public:

    int    order;
    double l1;
    double l2;
    double rotation = 0.0;
    boost::random::mt19937                     rng;
    boost::random::uniform_real_distribution<> phi_;
    boost::random::uniform_real_distribution<> psi_;
    Epetra_Map    * CellsMap;
    Epetra_Vector * GaussianRF;

    shinozuka_layeredcomp_2d(int & nu);
    shinozuka_layeredcomp_2d(int & nu, double & L1, double & L2);
    ~shinozuka_layeredcomp_2d();

    template<typename typearg>
    double tau_beta(typearg & beta);
    double s_tau(double & tau);

    void generator_gauss_points(Epetra_SerialDenseVector & v, mesh & Mesh, std::vector<int> & phase);
    void generator_one_gauss_point(Epetra_SerialDenseVector & v, mesh & Mesh, std::vector<int> & phase, double & xi, double & eta, double & zeta);

    void icdf_gamma(Epetra_Vector & V, Epetra_Vector & G, double & alpha, double & beta);
    void icdf_beta(Epetra_Vector & V, Epetra_Vector & B, double & tau1, double & tau2);
};
#endif
