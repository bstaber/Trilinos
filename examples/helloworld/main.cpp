#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

int main(int argc, char *argv[]){

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

    std::cout << "Hello from " << Comm.MyPID() << " out of " << Comm.NumProc() << "\n";
    Comm.Barrier();
    if (Comm.MyPID()==0){
	std::cout << "EXIT_SUCCESS\n";
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
