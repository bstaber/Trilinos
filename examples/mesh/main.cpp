#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "meshpp.hpp"

int main(int argc, char *argv[]){
	
#ifdef HAVE_MPI
MPI_Init(&argc, &argv);
Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
Epetra_SerialComm Comm;
#endif
    
    std::string mesh_file = "/Users/brian/Documents/GitHub/Trilinos/nrl/mesh/composite_hexa_32.msh";
    mesh Mesh(Comm, mesh_file);
    
    if (Comm.MyPID()==0){
        for (unsigned int i=0; i<Mesh.n_cells; ++i){
            std::cout << Mesh.epart[i] << "\n";
        }
    }
    
#ifdef HAVE_MPI
MPI_Finalize();
#endif
return 0;
    
}
