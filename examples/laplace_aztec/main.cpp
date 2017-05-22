#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "meshpp.hpp"
#include "laplacepp.hpp"

int main(int argc, char *argv[]){
	
#ifdef HAVE_MPI
MPI_Init(&argc, &argv);
Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
Epetra_SerialComm Comm;
#endif
    
    //Create a parameter list and two sublists "Mesh" and "Laplace".
    Teuchos::ParameterList Laplace_parameters;
    Teuchos::ParameterList & mesh_sublist = Laplace_parameters.sublist("Mesh");
    Teuchos::ParameterList & aztec_sublist = Laplace_parameters.sublist("Aztec");
    Teuchos::ParameterList & amesos_sublist = Laplace_parameters.sublist("Amesos");
    
    unsigned int nb_phys_groups = 4;
    mesh_sublist.set("mesh_file","test2.msh");
    mesh_sublist.set("boundary_file","nodes_to_boundaries_less_bumpy.txt");
    mesh_sublist.set("nb_phys_groups",nb_phys_groups);
    
    aztec_sublist.set("solver", "cg");
    aztec_sublist.set("precond", "jacobi");
    //aztec_sublist.set("subdomain_solve", "ilut");
    aztec_sublist.set("AZ_conv", "noscaled");
    aztec_sublist.set("AZ_tol", 1e-6);
    aztec_sublist.set("AZ_output", 0);
    aztec_sublist.set("AZ_diagnostics", 0);
    aztec_sublist.set("AZ_reorder", 1);
    
    bool display = 1;
    amesos_sublist.set("display",display);
    amesos_sublist.set("solver_type","Mumps");
    
    //$$ Laplace example
    Teuchos::RCP<laplace> LaplaceO = Teuchos::rcp(new laplace(Comm,Laplace_parameters.sublist("Mesh")));
    Epetra_FECrsMatrix matrix(Copy,*LaplaceO->FEGraph);
    Epetra_Vector lhs(*LaplaceO->StandardMap);
    Epetra_FEVector rhs(*LaplaceO->StandardMap);
    //Problem number one with aztec
    int bc_indx[2];
    double bc_val[2];
    bc_indx[0] = 0; bc_indx[1] = 1;
    bc_val[0] = 0.0; bc_val[1] = 1.0;
    LaplaceO->solve_aztec(Laplace_parameters.sublist("Aztec"), matrix, lhs, rhs, &bc_indx[0], &bc_val[0]);
    LaplaceO->print_solution(lhs, "laplace_inner_to_outer_aztec.mtx");
    //Problem number two with aztec
    bc_indx[0] = 2; bc_indx[1] = 3;
    LaplaceO->solve_aztec(Laplace_parameters.sublist("Aztec"), matrix, lhs, rhs, &bc_indx[0], &bc_val[0]);
    LaplaceO->print_solution(lhs, "laplace_inlet_to_outlet_aztec.mtx");
    //$$ End example
    
#ifdef HAVE_MPI
MPI_Finalize();
#endif
return 0;
    
}
