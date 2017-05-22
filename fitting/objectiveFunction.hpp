#define _USE_MATH_DEFINES
#include <math.h>

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

#include <stdio.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseSolver.h"
#include <random>

template<class Real, class XPrim=ROL::StdVector<Real>, class XDual=ROL::StdVector<Real> >
class objectiveFunction : public ROL::Objective<Real>{
    
//    typedef std::vector<Real> vector;
//    typedef Vector<Real>      V;
//    typedef typename vector::size_type uint;
    
private:
    
    template<class VectorType>
    Teuchos::RCP< const std::vector<Real> > getVector( const ROL::Vector<Real>& x ) {
        return Teuchos::dyn_cast<const VectorType>((x)).getVector();
    }

    template<class VectorType>
    Teuchos::RCP< std::vector<Real> > getVector( ROL::Vector<Real>& x ) {
        return Teuchos::dyn_cast<VectorType>(x).getVector();
    }
    
public:
    
    Epetra_SerialDenseVector y1_exp_circ, y1_exp_long;
    Epetra_SerialDenseVector tau11_exp_circ, tau11_exp_long;
    int ncirc, nlong;
    
    objectiveFunction(std::string filename_circ, std::string filename_long){
        set_experimental_data(filename_circ,filename_long);
    }
    objectiveFunction(){
    }
    //~objectiveFunction();
    
    Real value(const ROL::Vector<Real> &x, Real &tol){
        
        Teuchos::RCP<const std::vector<Real> > xp = getVector<XPrim>(x);
        uint n = xp->size();
        
        Epetra_SerialDenseVector p(n);
        for (unsigned int i=0; i<n; ++i){
            p(i) = (*xp)[i];
        }
        
        Real val = compute_error(p);
        return val;
    }
    
    double compute_error(Epetra_SerialDenseVector & p){
        
        double residual, tau;
        double ecirc = 0.0;
        double elong = 0.0;
        
        Epetra_SerialDenseVector y(3);
        
        y(1) = 1.0; y(2) = 1.0;
        for (unsigned int i=0; i<ncirc; ++i){
            y(0) = y1_exp_circ(i);
            residual = newton(y,p);
            compute_stress(y,p,tau);
            ecirc += (tau-tau11_exp_circ(i))*(tau-tau11_exp_circ(i));
            std::cout << std::setw(10) << y(0) << std::setw(10) <<  y(1) << std::setw(10) << y(2) << "\n";
        }
        ecirc/=(tau11_exp_circ.Norm2()*tau11_exp_circ.Norm2());
        
        p(6) = M_PI/2.0 - p(6);
        y(1) = 1.0; y(2) = 1.0;
        for (unsigned int j=0; j<nlong; ++j){
            y(0) = y1_exp_long(j);
            residual = newton(y,p);
            compute_stress(y,p,tau);
            elong += (tau-tau11_exp_long(j))*(tau-tau11_exp_long(j));
        }
        elong/=(tau11_exp_long.Norm2()*tau11_exp_long.Norm2());
        
        double kappa = 2.0*p(2)*p(4)*p(4);
        double mu = 2.0*p(0) + p(1)*3.0*std::sqrt(3.0);
        double nu = (3.0*kappa-2.0*mu)/(2.0*(3.0*kappa+mu));
        double enu = (nu-0.46)*(nu-0.46)/(0.46*0.46);
        
        double error = ecirc + elong + enu;
        return error;
    }
    
    double newton(Epetra_SerialDenseVector & y, Epetra_SerialDenseVector & xp){
        
        Epetra_SerialDenseSolver solver;
        
        double rhs_inf = 1.0;
        int iter = 0;
        
        Epetra_SerialDenseVector f(2);
        Epetra_SerialDenseVector dx(2);
        Epetra_SerialDenseMatrix A(2,2);
        
        while (rhs_inf>1e-10){
            
            compute_jacobian_and_residual(y,xp,f,A);
            
            iter++;
            if (iter>1){
                rhs_inf = f.NormInf();
                //std::cout << std::setw(15) << rhs_inf << std::setw(15) << y(0) << std::setw(15) << y(1) << std::setw(15) << y(2) << "\n";
                if (rhs_inf>1e8){
                    std::cout << "Newton-Raphson failed to converge.\n";
                    break;
                }
            }
            
            
            if (iter>1000){
                //std::cout << "Newton-Raphson: Iteration number exceeds 1000.";
                /*for (unsigned int i=0; i<xp.Length(); ++i){
                    std::cout << std::setw(10) << xp(i);
                }
                std::cout << "\n";*/
                break;
            }
            
            solver.SetMatrix(A);
            solver.SetVectors(dx,f);
            int error = solver.Solve();
            /*if (error){
             std::cout << "Inverse Isoparametric Mapping: Error " << error << " with Epetra_SerialDenseSolver at iteration " << iter << ".\n";
             }*/
            y(1) += dx(0);
            y(2) += dx(1);
        }

        return rhs_inf;
    }
    
    void compute_jacobian_and_residual(Epetra_SerialDenseVector & y, Epetra_SerialDenseVector & xp, Epetra_SerialDenseVector & f, Epetra_SerialDenseMatrix & A){
        
        double mu1 = xp(0); double mu2 = xp(1); double mu3 = xp(2); double mu4 = xp(3);
        double beta3 = xp(4); double beta4 = xp(5);
        double theta = xp(6); double c = cos(theta); double s = sin(theta);
        
        double J = y(0)*y(1)*y(2);
        double alpha = std::pow(J,-2.0/3.0);
        double beta = 1.0/(J*J);
        
        double I1 = y(0)*y(0) + y(1)*y(1) + y(2)*y(2);
        double I2 = y(0)*y(0)*y(1)*y(1) + y(0)*y(0)*y(2)*y(2) + y(1)*y(1)*y(2)*y(2);
        double I4 = y(0)*y(0)*c*c + y(1)*y(1)*s*s;
        double I5 = y(0)*y(0)*y(0)*y(0)*c*c + y(1)*y(1)*y(1)*y(1)*s*s;
        double K3 = I1*I4-I5;
        
        double pressure = mu3*beta3*( std::pow(J,beta3-1.0) - std::pow(J,-beta3-1.0) );
        double dpressure = mu3*beta3*( (beta3-1.0)*std::pow(J,beta3-2.0) + (beta3+1.0)*std::pow(J,-beta3-2.0) );
        
        f(0) = 2.0*mu1*alpha*(y(1)*y(1)-(1.0/3.0)*I1) + 2.0*mu2*beta*((3.0/2.0)*std::sqrt(I2)*(I1*y(1)*y(1)-y(1)*y(1)*y(1)*y(1)) - I2*std::sqrt(I2)) + J*pressure;
        f(1) = 2.0*mu1*alpha*(y(2)*y(2)-(1.0/3.0)*I1) + 2.0*mu2*beta*((3.0/2.0)*std::sqrt(I2)*(I1*y(2)*y(2)-y(2)*y(2)*y(2)*y(2)) - I2*std::sqrt(I2)) + J*pressure;
        
        A(0,0) = -(4.0/3.0)*(1.0/y(1))*mu1*alpha*(y(1)*y(1)-(1.0/3.0)*I1) + 2.0*mu1*alpha*(2.0*y(1)-(2.0/3.0)*y(1))
                 - 4.0*beta*mu2*(1.0/y(1))*( (3.0/2.0)*std::sqrt(I2)*(I1*y(1)*y(1)-y(1)*y(1)*y(1)*y(1)) - std::sqrt(I2)*I2 )
                 + 2.0*mu2*beta*( (3.0/2.0)*(1.0/std::sqrt(I2))*y(1)*(y(0)*y(0)+y(2)*y(2))*(I1*y(1)*y(1)-y(1)*y(1)*y(1)*y(1)) + (3.0/2.0)*std::sqrt(I2)*(-2.0*y(1)*y(1)*y(1)+2.0*y(1)*I1) - 3.0*std::sqrt(I2)*y(1)*(y(0)*y(0)+y(2)*y(2)) )
                 + J*(1.0/y(1))*(pressure + J*dpressure);
        
        A(0,1) = -(4.0/3.0)*(1.0/y(2))*mu1*alpha*(y(1)*y(1)-(1.0/3.0)*I1) - (4.0/3.0)*mu1*alpha*y(2)
                 - 4.0*mu2*beta*(1.0/y(2))*( (3.0/2.0)*std::sqrt(I2)*(I1*y(1)*y(1)-y(1)*y(1)*y(1)*y(1))-std::sqrt(I2)*I2 )
                 + 2.0*mu2*beta*( (3.0/2.0)*(1.0/std::sqrt(I2))*y(2)*(y(0)*y(0)+y(1)*y(1))*(I1*y(1)*y(1)-y(1)*y(1)*y(1)*y(1)) + 3.0*std::sqrt(I2)*y(1)*y(1)*y(2) - 3.0*std::sqrt(I2)*y(2)*(y(0)*y(0)+y(1)*y(1)) )
                 + J*(1.0/y(2))*(pressure + J*dpressure);
        
        A(1,0) = -(4.0/3.0)*mu1*alpha*(1.0/y(1))*(y(2)*y(2)-(1.0/3.0)*I1) - (4.0/3.0)*mu1*alpha*y(1)
                 - 4.0*mu2*beta*(1.0/y(1))*( (3.0/2.0)*std::sqrt(I2)*(I1*y(2)*y(2)-y(2)*y(2)*y(2)*y(2)) - std::sqrt(I2)*I2 )
                 + 2.0*mu2*beta*( (3.0/2.0)*(1.0/std::sqrt(I2))*y(1)*(y(0)*y(0)+y(2)*y(2))*(I1*y(2)*y(2)-y(2)*y(2)*y(2)*y(2)) + 3.0*std::sqrt(I2)*y(1)*y(2)*y(2) - 3.0*std::sqrt(I2)*y(1)*(y(0)*y(0)+y(2)*y(2)) )
                 + J*(1.0/y(1))*(pressure + J*dpressure);
        
        A(1,1) = -(4.0/3.0)*mu1*alpha*(1.0/y(2))*(y(2)*y(2)-(1.0/3.0)*I1) + 2.0*mu1*alpha*(2.0*y(2)-(2.0/3.0)*y(2))
                 - 4.0*mu2*beta*(1.0/y(2))*( (3.0/2.0)*std::sqrt(I2)*(I1*y(2)*y(2)-y(2)*y(2)*y(2)*y(2))-std::sqrt(I2)*I2 )
                 + 2.0*mu2*beta*( (3.0/2.0)*(1.0/std::sqrt(I2))*y(2)*(y(0)*y(0)+y(1)*y(1))*(I1*y(2)*y(2)-y(2)*y(2)*y(2)*y(2)) + (3.0/2.0)*std::sqrt(I2)*(-2.0*y(2)*y(2)*y(2)+2.0*y(2)*I1) - 3.0*std::sqrt(I2)*y(2)*(y(0)*y(0)+y(1)*y(1)) )
                 + J*(1.0/y(2))*(pressure + J*dpressure);
        
        
        if (I4>1.0){
            f(0) += 8.0*mu4*(I4-1.0)*exp(beta4*(I4-1.0)*(I4-1.0))*y(1)*y(1)*s*s;
            A(0,0) += 8.0*mu4*exp(beta4*(I4-1.0)*(I4-1.0))*y(1)*y(1)*s*s*y(1)*s + 16.0*mu4*beta4*(I4-1.0)*(I4-1.0)*exp(beta4*(I4-1.0)*(I4-1.0))*y(1)*y(1)*y(1)*y(1)*s*s*s*s + 16.0*mu4*(I4-1.0)*(I4-1.0)*exp(beta4*(I4-1.0)*(I4-1.0))*y(1)*s*s;
        }
         
        /*if (K3>2.0){
            f(0) += 4.0*mu4*beta4*std::pow(K3-2.0,beta4-1.0)*(I4 + s*s*(I1-2.0*y(1)*y(1)))*y(1)*y(1);
            f(1) += 4.0*mu4*beta4*std::pow(K3-2.0,beta4-1.0)*I4*y(2)*y(2);
            
            A(0,0) += 4.0*mu4*beta4*(beta4-1.0)*std::pow(K3-2.0,beta4-2.0)*(2.0*y(1)*I4+2.0*I1*y(1)*s*s-4.0*y(1)*y(1)*y(1)*s*s)*(I4+s*s*(I1-2.0*y(1)*y(1)))*y(1)*y(1) + 4.0*mu4*beta4*std::pow(K3-2.0,beta4-1.0)*(2.0*y(1)*s*s+s*s*(2.0*y(1)-4.0*y(1)))*y(1)*y(1) + 8.0*mu4*beta4*std::pow(K3-2.0,beta4-1.0)*(I4+s*s*(I1-2.0*y(1)*y(1)))*y(1);
            A(0,1) += 4.0*mu4*beta4*(beta4-1.0)*std::pow(K3-2.0,beta4-2.0)*(I4+s*s*(I1-2.0*y(1)*y(1)))*y(1)*y(1)*2.0*y(2)*I4 + 8.0*mu4*beta4*std::pow(K3-2.0,beta4-1.0)*y(1)*y(1)*y(2)*s*s;
            A(1,0) += 4.0*mu4*beta4*(beta4-1.0)*std::pow(K3-2.0,beta4-2.0)*I4*y(2)*y(2)*(2.0*y(1)*I4+2.0*I1*y(1)*s*s-4.0*y(1)*y(1)*y(1)*s*s) + 8.0*mu4*beta4*std::pow(K3-2.0,beta4-1.0)*I4*y(2);
            A(1,1) += 4.0*mu4*beta4*(beta4-1.0)*std::pow(K3-2.0,beta4-2.0)*I4*y(2)*y(2)*2.0*y(2)*I4 + 8.0*mu4*beta4*std::pow(K3-2.0,beta4-1.0)*I4*y(2);
        }*/
        
        f.Scale(-1.0);
    }
    
    void compute_stress(Epetra_SerialDenseVector & y, Epetra_SerialDenseVector & xp, double & tau){
        
        double mu1 = xp(0); double mu2 = xp(1); double mu3 = xp(2); double mu4 = xp(3);
        double beta3 = xp(4); double beta4 = xp(5);
        double theta = xp(6); double c = cos(theta); double s = sin(theta);
        
        double J = y(0)*y(1)*y(2);
        double alpha = std::pow(J,-2.0/3.0);
        double beta = 1.0/(J*J);
        
        double I1 = y(0)*y(0) + y(1)*y(1) + y(2)*y(2);
        double I2 = y(0)*y(0)*y(1)*y(1) + y(0)*y(0)*y(2)*y(2) + y(1)*y(1)*y(2)*y(2);
        double I4 = y(0)*y(0)*c*c + y(1)*y(1)*s*s;
        double I5 = y(0)*y(0)*y(0)*y(0)*c*c + y(1)*y(1)*y(1)*y(1)*s*s;
        double K3 = I1*I4-I5;
        
        double pressure = mu3*beta3*( std::pow(J,beta3-1.0) - std::pow(J,-beta3-1.0) );
        
        tau = 2.0*mu1*alpha*(y(0)*y(0)-(1.0/3.0)*I1) + 2.0*mu2*beta*((3.0/2.0)*std::sqrt(I2)*(I1*y(0)*y(0)-y(0)*y(0)*y(0)*y(0))-std::sqrt(I2)*I2) + J*pressure;
        if (K3>2.0){
            tau += 8.0*mu4*(I4-1.0)*exp(beta4*(I4-1.0)*(I4-1.0))*y(0)*y(0)*c*c;
            //tau += 4.0*mu4*beta4*std::pow(K3-2.0,beta4-1.0)*(I4+c*c*(I1-2.0*y(0)*y(0)))*y(0)*y(0);
        }
        
    }
    
    void set_experimental_data(std::string filename_circ, std::string filename_long){
        
        double stretch, stress;
        int ndata_circ, ndata_long;
        std::ifstream file1, file2;
        
        std::cout << std::setw(30) << "BEGIN: READ DATA BASE\n";
        std::cout << "--------------------------------------------\n";
        file1.open(filename_circ);
        file2.open(filename_long);
        if (file1.is_open() && file2.is_open()){
            file1 >> ndata_circ;
            y1_exp_circ.Resize(ndata_circ);
            tau11_exp_circ.Resize(ndata_circ);
            for (unsigned int i=0; i<ndata_circ; ++i){
                file1 >> stretch;
                file1 >> stress;
                y1_exp_circ(i) = stretch;
                tau11_exp_circ(i) = stress;
            }
            file2 >> ndata_long;
            y1_exp_long.Resize(ndata_long);
            tau11_exp_long.Resize(ndata_long);
            for (unsigned int i=0; i<ndata_long; ++i){
                file2 >> stretch;
                file2 >> stress;
                y1_exp_long(i) = stretch;
                tau11_exp_long(i) = stress;
            }
            file1.close();
            file2.close();
        }
        else{
            std::cout << "Couldn't open something\n";
        }
        ncirc = y1_exp_circ.Length();
        nlong = y1_exp_long.Length();
        
        std::cout << std::setw(10) << "y1" << std::setw(10) << "tau11" << std::setw(10) << "y1" << std::setw(10) << "tau11\n";
        for (unsigned int i=0; i<fmin(ncirc,nlong); ++i){
            std::cout << std::setw(10) << y1_exp_circ(i) << std::setw(10) << tau11_exp_circ(i) << std::setw(10) << y1_exp_long(i) << std::setw(10) << tau11_exp_long(i) << "\n";
        }
        std::cout << std::setw(30) << "END: READ DATA BASE\n";
        std::cout << "--------------------------------------------\n";
        
    }
    
};
