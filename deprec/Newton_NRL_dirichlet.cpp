#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "NRL_ModelF.hpp"
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
    aztec.set("solver", "gmres");
    aztec.set("kspace", 500);
    aztec.set("orthog", "classical");
    //aztec.set("precond", "jacobi");
    aztec.set("precond", "dom_decomp");
    aztec.set("subdomain_solve", "ilut");
    aztec.set("overlap", 2);
    aztec.set("type_overlap", "symmetric");
    aztec.set("pre_calc", "recalc");
    aztec.set("ilut_fill",3.0);
    aztec.set("athresh", 0.0);
    aztec.set("rthresh", 0.0);
    aztec.set("drop",0.0);
    aztec.set("AZ_conv", "noscaled");
    aztec.set("AZ_tol", 1e-6);
    aztec.set("AZ_output", 0);
    aztec.set("AZ_diagnostics", "all");
    aztec.set("AZ_reorder", 1);
    
    Teuchos::ParameterList & newton = Parameters.sublist("Newton");
    
    newton.set("delta", 0.05);
    newton.set("iterMax", 10);
    newton.set("nbBisMax", 5);
    newton.set("NormFTol", 1e-6);
    newton.set("NormFMax", 1e7);
    newton.set("eps", 1e-8);
    newton.set("success_parameter", 1.0);
    newton.set("failure_parameter", 2.0);
    newton.set("number_of_loads", 1);
    newton.set("bc_disp",5.0);
    newton.set("pressure_load",0.0);

    Teuchos::ParameterList & mesh_sublist = Parameters.sublist("Mesh");
    //mesh_sublist.set("mesh_file","/Users/Brian/Documents/Thesis/0-Trilinos/Meshes/NRL_composite/composite.msh");
    mesh_sublist.set("mesh_file","/home/s/staber/Trilinos/NRL_composite/composite_hexa.msh");
    mesh_sublist.set("n_ply",8);
    //unsigned int nb_phys_groups = 4;
    //mesh_sublist.set("boundary_file","nodes_to_boundaries_less_bumpy.txt");
    //mesh_sublist.set("nb_phys_groups",nb_phys_groups);
    
    Teuchos::ParameterList & model_B = Parameters.sublist("Model_B");
    Teuchos::ParameterList & model_F = Parameters.sublist("Model_F");
    //Model B
    model_B.set("mu1",8.0);
    model_B.set("mu2",0.0);
    model_B.set("epsilon1",10.0);
    model_B.set("mu4",3000.0);
    model_B.set("beta4",3.0);
    model_B.set("beta5",2.0);
    model_B.set("angle",0.26);
    //Model F
    model_F.set("m1",5.357);
    model_F.set("m2",7.923);
    model_F.set("beta3",-1.0/2.0);
    model_F.set("beta4",2.495);
    model_F.set("beta5",0.694);
    model_F.set("angle",0.26);
    
    Teuchos::RCP<Interface_dirichlet> my_interface = Teuchos::rcp(new Interface_dirichlet(Comm,Parameters));
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*my_interface,Parameters));
    
    Newton->Initialization();
    
    int error1 = Newton->Krylov_Method();
    double xi = 1.0/2.0;
    my_interface->compute_green_lagrange(*Newton->x,xi,xi,xi);
    
#ifdef HAVE_MPI
MPI_Finalize();
#endif
return 0;
    
}
