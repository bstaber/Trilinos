#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//#include "Arteries_ModelA.hpp"
//#include "Arteries_ModelB.hpp"
//#include "Arteries_ModelC.hpp"
//#include "Arteries_ModelC.hpp"
//#include "Arteries_ModelC_slipbc.hpp"
//#include "Arteries_ModelC_slipbc_and_clamp.hpp"
//#include "Arteries_ModelC_semislipbc.hpp"
//#include "Arteries_ModelC_deterministic_neumann.hpp"
#include "Arteries_ModelC_deterministic_neumann.hpp"
#include "Newton_Raphsonpp.hpp"

int main(int argc, char *argv[]){
    
#ifdef HAVE_MPI
MPI_Init(&argc, &argv);
Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
Epetra_SerialComm Comm;
#endif
 
    Teuchos::ParameterList Parameters;
    Teuchos::ParameterList & aztec = Parameters.sublist("Krylov");
    Teuchos::ParameterList & newton = Parameters.sublist("Newton");
    Teuchos::ParameterList & mesh_sublist = Parameters.sublist("Mesh");
    Teuchos::ParameterList & laplace = Parameters.sublist("Laplace");
    Teuchos::ParameterList & modelc = Parameters.sublist("ModelC");
    
    aztec.set("solver","gmres");
    aztec.set("kspace",500);
    aztec.set("orthog","classical");
    aztec.set("precond","dom_decomp");
    aztec.set("subdomain_solve","ilut");
    aztec.set("overlap",2);
    aztec.set("type_overlap","symmetric");
    aztec.set("pre_calc","recalc");
    aztec.set("ilut_fill",4.0);
    aztec.set("athresh",0.0);
    aztec.set("rthresh",0.0);
    aztec.set("drop",0.0);
    aztec.set("AZ_conv","noscaled");
    aztec.set("AZ_tol",1e-6);
    aztec.set("AZ_output",0);
    aztec.set("AZ_diagnostics","all");
    aztec.set("AZ_reorder",1);
    
    newton.set("delta",0.0025);
    newton.set("iterMax",10);
    newton.set("nbBisMax",5);
    newton.set("NormFTol",1e-6);
    newton.set("NormFMax",1e7);
    newton.set("eps",1e-8);
    newton.set("success_parameter",2.0);
    newton.set("failure_parameter",2.0);
    newton.set("number_of_loads",1);
    newton.set("bc_disp",0.0);
    newton.set("pressure_load",12.0*1000.0);

    unsigned int nb_phys_groups = 4;
    mesh_sublist.set("mesh_file","/home/s/staber/Trilinos/Meshes/Arterial_wall/media_flatboundaries.msh");
    mesh_sublist.set("boundary_file","/home/s/staber/Trilinos/Meshes/Arterial_wall/nodes_to_boundaries_media.txt");
    mesh_sublist.set("nb_phys_groups",nb_phys_groups);
    
    laplace.set("solver", "cg");
    laplace.set("precond", "dom_decomp");
    laplace.set("subdomain_solve", "icc");
    laplace.set("overlap",0);
    laplace.set("AZ_conv", "noscaled");
    laplace.set("AZ_tol", 1e-6);
    laplace.set("AZ_output", 0);
    laplace.set("AZ_diagnostics", 0);
    laplace.set("AZ_reorder", 1);
    
    modelc.set("mu1", 3.0759*1000.0);
    modelc.set("mu2", 0.21828*1000.0);
    modelc.set("mu3", 5.2622*1000.0);
    modelc.set("mu4", 37.259*1000.0);
    modelc.set("beta3", 3.5384);
    modelc.set("beta4", 496.49);
    modelc.set("theta", 0.82323);
        
    Teuchos::RCP<Interface_arteries> my_interface = Teuchos::rcp(new Interface_arteries(Comm,Parameters));
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*my_interface,Parameters));
    Newton->x->PutScalar(0.0);
    int error1 = Newton->Solve_with_Aztec();
    std::string name = "Newton_solution_slipbc.mtx";
    Newton->print_newton_solution(*Newton->x,name);
    
#ifdef HAVE_MPI
MPI_Finalize();
#endif
return 0;
    
}
