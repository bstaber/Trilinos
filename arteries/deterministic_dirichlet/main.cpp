#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Arteries_ModelC_deterministic_dirichlet.hpp"
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
    
    newton.set("delta",0.25);
    newton.set("iterMax",10);
    newton.set("nbBisMax",5);
    newton.set("NormFTol",1e-6);
    newton.set("NormFMax",1e7);
    newton.set("eps",1e-8);
    newton.set("success_parameter",2.0);
    newton.set("failure_parameter",2.0);
    newton.set("number_of_loads",1);
    newton.set("bc_disp",1.0);
    newton.set("pressure_load",0.0/1000.0);
    
    unsigned int nb_phys_groups = 4;
    //mesh_sublist.set("mesh_file","/Users/Brian/Documents/Thesis/0-Trilinos/Meshes/Arterial_wall/arterial_wall.msh");
    //mesh_sublist.set("boundary_file","/Users/Brian/Documents/Thesis/0-Trilinos/Meshes/Arterial_wall/nodes_to_boundaries_less_bumpy.txt");
    mesh_sublist.set("mesh_file","/home/s/staber/Trilinos/Meshes/Arterial_wall/media.msh"); //arterial_wall.msh");
    mesh_sublist.set("boundary_file","/home/s/staber/Trilinos/Meshes/Arterial_wall/nodes_to_boundaries_media.txt"); //nodes_to_boundaries.txt");
    mesh_sublist.set("nb_phys_groups",nb_phys_groups);
    
    Teuchos::ParameterList & laplace = Parameters.sublist("Laplace");
    laplace.set("solver", "cg");
    laplace.set("precond", "dom_decomp");
    laplace.set("subdomain_solve", "icc");
    laplace.set("overlap",0);
    laplace.set("AZ_conv", "noscaled");
    laplace.set("AZ_tol", 1e-6);
    laplace.set("AZ_output", 0);
    laplace.set("AZ_diagnostics", 0);
    laplace.set("AZ_reorder", 1);
    
    Teuchos::RCP<DirichletInletOutlet_PolyconvexHGO> my_interface = Teuchos::rcp(new DirichletInletOutlet_PolyconvexHGO(Comm,Parameters));
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*my_interface,Parameters));
    
    my_interface->fixed = 2;
    my_interface->moved = 3;
    my_interface->xxmax = 11.696;
    
    Newton->Initialization();
    
    int throwint;
    throwint = Newton->Solve_with_Aztec(true);
    
    my_interface->fixed = 3;
    my_interface->moved = 2;
    my_interface->xxmax = 0.15685;
    
    throwint = Newton->Solve_with_Aztec(true);
    
    std::string name="dirichlet_arteries.mtx";
    Newton->print_newton_solution(name);
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
    
}
