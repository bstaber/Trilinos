#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//Media + adventitia GMRF

//#include "Arteries_ModelA.hpp"
//#include "Arteries_ModelB.hpp"
//#include "Arteries_ModelC.hpp"
#include "Arteries_ModelC_GMRF.hpp"
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
    //mesh_sublist.set("mesh_file","/Users/Brian/Documents/Thesis/0-Trilinos/Projects/Arterial_wall/Mesh/arterial_wall.msh");
    //mesh_sublist.set("boundary_file","/Users/Brian/Documents/Thesis/0-Trilinos/Projects/Arterial_wall/Mesh/nodes_to_boundaries.txt");
    mesh_sublist.set("mesh_file",path+"arterial_wall.msh");
    mesh_sublist.set("boundary_file",path+"nodes_to_boundaries.txt");
    mesh_sublist.set("nb_phys_groups",nb_phys_groups);
    
    Teuchos::ParameterList & laplace = Parameters.sublist("Laplace");
    laplace.set("solver", "cg");
    laplace.set("precond", "jacobi");
    laplace.set("AZ_conv", "noscaled");
    laplace.set("AZ_tol", 1e-6);
    laplace.set("AZ_output", 0);
    laplace.set("AZ_diagnostics", 0);
    laplace.set("AZ_reorder", 1);
    
    Teuchos::RCP<Interface_arteries> my_interface = Teuchos::rcp(new Interface_arteries(Comm,Parameters));
    
    std::ifstream parameters_file_1, parameters_file_2, parameters_file_3, parameters_file_4, parameters_file_5, parameters_file_6, parameters_file_7, parameters_file_8;
    
    parameters_file_1.open(path+"C1_med.txt");
    parameters_file_2.open(path+"C2_med.txt");
    parameters_file_3.open(path+"U1_med.txt");
    parameters_file_4.open(path+"G4_med.txt");
    parameters_file_5.open(path+"C1_adv.txt");
    parameters_file_6.open(path+"C2_adv.txt");
    parameters_file_7.open(path+"U1_adv.txt");
    parameters_file_8.open(path+"G4_adv.txt");
    
    unsigned int n_cells_med = 140556;
    unsigned int n_nodes_med = 29600;
    unsigned int n_cells_adv = 70278;
    unsigned int n_nodes_adv = 17760;

    int error;
    my_interface->get_media_adventitia(n_cells_med,n_nodes_med,n_cells_adv,n_nodes_adv,path);
    if (parameters_file_1.is_open() && parameters_file_2.is_open() && parameters_file_3.is_open() && parameters_file_4.is_open() && parameters_file_5.is_open() && parameters_file_6.is_open() && parameters_file_7.is_open() && parameters_file_8.is_open()){
        
        for (unsigned nmc=0; nmc<6; ++nmc){
            for (int i=0; i<n_nodes_med; ++i){
                parameters_file_1 >> my_interface->c1_med(i);
                parameters_file_2 >> my_interface->c2_med(i);
                parameters_file_3 >> my_interface->u1_med(i);
                parameters_file_4 >> my_interface->mu4_med(i);
            }
            for (int j=0; j<n_nodes_adv; ++j){
                parameters_file_5 >> my_interface->c1_adv(j);
                parameters_file_6 >> my_interface->c2_adv(j);
                parameters_file_7 >> my_interface->u1_adv(j);
                parameters_file_8 >> my_interface->mu4_adv(j);
            }
            
            Teuchos::RCP<Newton_Raphson> Newton = Teuchos::rcp(new Newton_Raphson(*my_interface,Parameters));
            Newton->x->PutScalar(0.0);
            error = Newton->Krylov_Method();
        }
        
        Comm.Barrier();
        parameters_file_1.close();
        parameters_file_2.close();
        parameters_file_3.close();
        parameters_file_4.close();
        parameters_file_5.close();
        parameters_file_6.close();
        parameters_file_7.close();
        parameters_file_8.close();
    }
    else{
        std::cout << "Couldn't open one of the parameters_file.\n";
    }
    
#ifdef HAVE_MPI
MPI_Finalize();
#endif
return 0;
    
}
