#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "objectiveFunction.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_StatusTest.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include <random>

typedef double RealT;

int main(int argc, char *argv[]){
    
    int dim = 7;
    RealT alpha = 1.e-2;
    
    std::string    xmlInFileName = "";
    
    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setDocString("TO DO.");
    
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;
        return parse_return;
    }
    
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    if(xmlInFileName.length()) {
        Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(*parlist));
    }
    parlist->sublist("General").sublist("Krylov").set("Absolute Tolerance",1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Relative Tolerance",1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit",50);
    // PDAS parameters.
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-8);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-6);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit", 1);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",(alpha>0.0)?alpha:1.e-4);
    // Status test parameters.
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-12);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-14);
    parlist->sublist("Status Test").set("Iteration Limit",10000);
    
    std::string path = "/Users/Brian/Documents/Thesis/0-Trilinos/Trilinos/Fitting/data/";
    //std::string path = "/home/s/staber/Fitting/data/";
    
    for (unsigned int s=0; s<1; ++s){
        
        std::string circu = path + "specimen" + std::to_string(s) + "_circ.txt";
        std::string longi = path + "specimen" + std::to_string(s) + "_long.txt";
        
        Teuchos::RCP<objectiveFunction<double>> obj = Teuchos::rcp(new objectiveFunction<double>(circu,longi));
        
        Teuchos::RCP<std::vector<RealT>> l_rcp = Teuchos::rcp( new std::vector<RealT>(dim) );
        Teuchos::RCP<std::vector<RealT>> u_rcp = Teuchos::rcp( new std::vector<RealT>(dim) );
        
        Teuchos::RCP<ROL::Vector<RealT>> lo = Teuchos::rcp( new ROL::StdVector<RealT>(l_rcp) );
        Teuchos::RCP<ROL::Vector<RealT>> up = Teuchos::rcp( new ROL::StdVector<RealT>(u_rcp) );
        
        (*l_rcp)[0] = 0.01; (*l_rcp)[1] = 0.01; (*l_rcp)[2] = 0.01; (*l_rcp)[3] = 0.01;
        (*l_rcp)[4] = 2.0+1e-6; (*l_rcp)[5] = 2.0+1e-6; (*l_rcp)[6] = M_PI/6.0;
        
        (*u_rcp)[0] = 1000.0; (*u_rcp)[1] = 1000.0; (*u_rcp)[2] = 100.0; (*u_rcp)[3] = 50000.0;
        (*u_rcp)[4] = 100.0; (*u_rcp)[5] = 10.0; (*u_rcp)[6] = M_PI/3.0;
        
        Teuchos::RCP<std::vector<RealT> > x_rcp = Teuchos::rcp( new std::vector<RealT> (dim, 0.0) );
        /*(*x_rcp)[0] = 0.0051967; (*x_rcp)[1] = 0.55309; (*x_rcp)[2] = 7.1776; (*x_rcp)[3] = 48987.0;
        (*x_rcp)[4] = 9.7874; (*x_rcp)[5] = 6.0099; (*x_rcp)[6] = 0.74212;*/
        ROL::StdVector<RealT> x(x_rcp);
        
        std::mt19937 G(Comm.MyPID());
        std::uniform_real_distribution<double> mu1(1e-3,10.0-1e-3);
        std::uniform_real_distribution<double> mu2(1e-3,10.0-1e-3);
        std::uniform_real_distribution<double> mu3(1e-3,10.0-1e-3);
        std::uniform_real_distribution<double> mu4(1e-3,10.0-1e-3);
        std::uniform_real_distribution<double> beta3(2.0+1e-3,10.0-1e-3);
        std::uniform_real_distribution<double> beta4(2.0+1e-3,10.0-1e-3);
        std::uniform_real_distribution<double> theta(M_PI/6.0,M_PI/3.0);
        
        ROL::BoundConstraint<RealT> icon(lo,up);
        
        Epetra_SerialDenseVector p(dim);
        double fval;
        
        std::cout << std::setw(5) << "nmc" << std::setw(10) << "mu1" << std::setw(10) << "mu2" << std::setw(10) << "mu3" << std::setw(10) << "mu4" << std::setw(10) << "beta3" << std::setw(10) << "beta4" << std::setw(10) << "theta" << std::setw(10) << "residual" << "\n";
        
        for (unsigned int nmc=0; nmc<1000; ++nmc){
            
            (*x_rcp)[0] = mu1(G); (*x_rcp)[1] = mu2(G); (*x_rcp)[2] = mu2(G); (*x_rcp)[3] = mu4(G);
            (*x_rcp)[4] = beta3(G); (*x_rcp)[5] = beta4(G); (*x_rcp)[6] = theta(G);
            
            Teuchos::RCP<ROL::Algorithm<RealT> > algo = Teuchos::rcp(new ROL::Algorithm<RealT>("Trust Region",*parlist,false));
            algo->run(x, *obj, icon, true, std::cout);
            
            for (unsigned int j=0; j<dim; ++j){
                p(j) = (*x_rcp)[j];
            }
            fval = obj->compute_error(p);
            
            std::cout << std::setw(5) << nmc;
            for (unsigned int j=0; j<dim; ++j){
                std::cout << std::setw(10) << p(j);
            }
            std::cout << std::setw(10) << fval;
            std::cout << "\n";
        }
        
    }
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
