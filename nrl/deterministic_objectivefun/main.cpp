#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Newton_Raphsonpp.hpp"
#include "objectiveFunction.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

int main(int argc, char *argv[]){
    
    std::string    xmlInFileName = "";
    std::string    xmlExtraFileName = "";
    
    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setOption("extra-file",&xmlExtraFileName,"Extra XML file to read into a parameter list");
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
    
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(new Teuchos::ParameterList);
    if(xmlInFileName.length()){
        Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(*paramList));
    }
    if(xmlExtraFileName.length()){
        Teuchos::updateParametersFromXmlFile(xmlExtraFileName, inoutArg(*paramList));
    }
    if (Comm.MyPID()==0){
        paramList->print(std::cout,2,true,true);
    }
    
    Teuchos::RCP<objectiveFunction<double>> obj = Teuchos::rcp(new objectiveFunction<double>(Comm,*paramList));
    
    /*Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    if(xmlInFileName.length()) {
        Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(*parlist));
    }
    parlist->sublist("General").sublist("Krylov").set("Absolute Tolerance",1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Relative Tolerance",1.e-4);
    parlist->sublist("General").sublist("Krylov").set("Iteration Limit",50);
    // PDAS parameters.
    double alpha = 1.e-2;
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Step Tolerance",1.e-8);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Relative Gradient Tolerance",1.e-6);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Iteration Limit", 1);
    parlist->sublist("Step").sublist("Primal Dual Active Set").set("Dual Scaling",(alpha>0.0)?alpha:1.e-4);
    // Status test parameters.
    parlist->sublist("Status Test").set("Gradient Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Step Tolerance",1.e-8);
    parlist->sublist("Status Test").set("Iteration Limit",10000);
    
    bool printHeader;
    if (Comm.MyPID()==0){
        printHeader = true;
    }
    else{
        printHeader = false;
    }
    Teuchos::RCP<std::vector<double>> l_rcp = Teuchos::rcp( new std::vector<double>(5) );
    Teuchos::RCP<std::vector<double>> u_rcp = Teuchos::rcp( new std::vector<double>(5) );
    
    Teuchos::RCP<ROL::Vector<double>> lo = Teuchos::rcp( new ROL::StdVector<double>(l_rcp) );
    Teuchos::RCP<ROL::Vector<double>> up = Teuchos::rcp( new ROL::StdVector<double>(u_rcp) );
    
    (*l_rcp)[0] = 1.e1; (*l_rcp)[1] = 1.e1; (*l_rcp)[2] = -1.0/2.0; (*l_rcp)[3] = 1.e-3; (*l_rcp)[4] = 1.e-3;
    (*u_rcp)[0] = 1.e5; (*u_rcp)[1] = 1.e5; (*u_rcp)[2] = 1.e1; (*u_rcp)[3] = 1.e1; (*u_rcp)[4] = 1.e1;
    
    ROL::BoundConstraint<double> icon(lo,up);
    
    Comm.Barrier();
    if (Comm.MyPID()==0){
        std::cout << "\n";
        std::cout << std::setw(5) << "nmc" << std::setw(20) << "value" << std::setw(20) << "m1" << std::setw(20) << "m2" << std::setw(20) << "beta3" << std::setw(20) << "beta4" << std::setw(20) << "beta5" << "\n";
    }
    
    boost::random::mt19937 rng;
    for (unsigned int nmc=0; nmc<1000; ++nmc){
        
        Teuchos::RCP<ROL::Algorithm<double> > algo =
        Teuchos::rcp(new ROL::Algorithm<double>("Trust Region",*parlist,false));
        
        Teuchos::RCP<std::vector<double> > x_rcp = Teuchos::rcp( new std::vector<double> (5, 0.0) );
     //(*x_rcp)[0] = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"m1");
     //(*x_rcp)[1] = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"m2");
     //(*x_rcp)[2] = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"beta3");
     //(*x_rcp)[3] = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"beta4");
     //(*x_rcp)[4] = Teuchos::getParameter<double>(paramList->sublist("ModelF"),"beta5");
    
        boost::random::uniform_real_distribution<> m1((*l_rcp)[0],(*u_rcp)[0]);
        boost::random::uniform_real_distribution<> m2((*l_rcp)[1],(*u_rcp)[1]);
        boost::random::uniform_real_distribution<> beta3((*l_rcp)[2],(*u_rcp)[2]);
        boost::random::uniform_real_distribution<> beta4((*l_rcp)[3],(*u_rcp)[3]);
        boost::random::uniform_real_distribution<> beta5((*l_rcp)[4],(*u_rcp)[4]);
        (*x_rcp)[0] = m1(rng);
        (*x_rcp)[1] = m2(rng);
        (*x_rcp)[2] = beta3(rng);
        (*x_rcp)[3] = beta4(rng);
        (*x_rcp)[4] = beta5(rng);
        ROL::StdVector<double> x(x_rcp);

        //algo->run(x, *obj, icon, printHeader, std::cout);
        double tol = 1e-6;
        double value = obj->value(x,tol);
    
        if (Comm.MyPID()==0){
            std::cout << std::setw(5) << nmc << std::setw(20) << value;
            for (unsigned int j=0; j<5; ++j){
                std::cout << std::setw(20) << (*x_rcp)[j];
            }
            std::cout << "\n";
        }
    }*/
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
