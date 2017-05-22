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
#include "Arteries_ModelC_GMRF_media.hpp"
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
    
    aztec.set("solver", "gmres");
    aztec.set("kspace", 500);
    aztec.set("orthog", "classical");
    aztec.set("precond", "dom_decomp");
    aztec.set("subdomain_solve", "ilut");
    aztec.set("overlap", 2);
    aztec.set("type_overlap", "symmetric");
    aztec.set("pre_calc", "recalc");
    aztec.set("ilut_fill", 4.0);
    aztec.set("athresh", 0.0);
    aztec.set("rthresh", 0.0);
    aztec.set("drop",0.0);
    aztec.set("AZ_conv", "noscaled");
    aztec.set("AZ_tol", 1e-6);
    aztec.set("AZ_output", 0);
    aztec.set("AZ_diagnostics", "all");
    aztec.set("AZ_reorder", 1);
    
    newton.set("delta", 0.005);
    newton.set("iterMax", 10);
    newton.set("nbBisMax", 5);
    newton.set("NormFTol", 1e-6);
    newton.set("NormFMax", 1e7);
    newton.set("eps", 1e-8);
    newton.set("success_parameter", 2.0);
    newton.set("failure_parameter", 2.0);
    newton.set("number_of_loads", 1);
    newton.set("bc_disp",0.0);
    newton.set("pressure_load",12.0/1000.0);

    std::string path = "/home/s/staber/Trilinos/Meshes/Arterial_wall/";
    unsigned int nb_phys_groups = 4;
    mesh_sublist.set("mesh_file",path+"media.msh");
    mesh_sublist.set("boundary_file",path+"nodes_to_boundaries_media.txt");
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
    
    Teuchos::RCP<Interface_arteries> my_interface = Teuchos::rcp(new Interface_arteries(Comm,Parameters));
    Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*my_interface,Parameters));
    
    std::ifstream parameters_file_1, parameters_file_2, parameters_file_3, parameters_file_4;
    
    parameters_file_1.open(path+"C1.txt");
    parameters_file_2.open(path+"C2.txt");
    parameters_file_3.open(path+"U1.txt");
    parameters_file_4.open(path+"G4.txt");
    
    unsigned int n_cells_p1_med = 297828;
    unsigned int n_nodes_p1_med = 58464;

    int error;
    my_interface->get_media(n_cells_p1_med,n_nodes_p1_med,path);
    if (parameters_file_1.is_open() && parameters_file_2.is_open() && parameters_file_3.is_open() && parameters_file_4.is_open()){
        
        for (unsigned nmc=0; nmc<6; ++nmc){
            for (int i=0; i<n_nodes_p1_med; ++i){
                parameters_file_1 >> my_interface->c1_gmrf(i);
                parameters_file_2 >> my_interface->c2_gmrf(i);
                parameters_file_3 >> my_interface->u1_gmrf(i);
                parameters_file_4 >> my_interface->mu4_gmrf(i);
            }
            Newton->x->PutScalar(0.0);
            error = Newton->Krylov_Method();
            std::string name="realization" + std::to_string(nmc) + ".mtx";
            Newton->print_newton_solution(*Newton->x,name);
        }
        Comm.Barrier();
        parameters_file_1.close();
        parameters_file_2.close();
        parameters_file_3.close();
        parameters_file_4.close();
    }
    else{
        std::cout << "Couldn't open one of the parameters_file.\n";
    }
    
#ifdef HAVE_MPI
MPI_Finalize();
#endif
return 0;
    
}
