#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "readnrldata.hpp"
#include "distributenrldata.hpp"

int main(int argc, char *argv[]){

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

    //std::string filemesh = "/Users/brian/Documents/GitHub/Trilinos/nrl/mesh/composite_hexa_32.msh";
    std::string filemesh = "/home/s/staber/Trilinos/nrl/mesh/composite_hexa_32.msh";
    mesh Mesh(Comm,filemesh, 1.0);
    std::string pathnrl = "/home/s/staber/Trilinos_results/nrl/data/";
    Teuchos::RCP<distributenrldata> nrldata = Teuchos::rcp(new distributenrldata(Mesh,pathnrl));

    if (Comm.MyPID()==0){
        std::cout << "MyPID" << std::setw(10) << "ncells" << std::setw(10) << "npoints" << "\n";
    }
    Comm.Barrier();
    std::cout << Comm.MyPID() << std::setw(10) << nrldata->local_cells.size() << std::setw(10) << nrldata->local_xi.size() << "\n";
    Comm.Barrier();
    /*if (Comm.MyPID()==1){
        for (unsigned int i=0; i<nrldata->local_cells.size(); ++i){
            std::cout << nrldata->local_xi[i] << std::setw(10) << nrldata->local_eta[i] << "\n";
        }
    }
    Comm.Barrier();
    */
    if (Comm.MyPID()==0){
        std::cout << nrldata->boundaryconditions << "\n";
        std::cout << nrldata->energy << "\n";
    }

    //Teuchos::RCP<readnrldata> nrldata = Teuchos::rcp(new readnrldata());
    /*if (Comm.MyPID()==0){
        std::cout << "npoints = " << nrldata->npoints << std::setw(10) << "nloads = " << nrldata->nloads << "\n";

        for (unsigned int i=0; i<nrldata->points.M(); ++i){
            std::cout << std::setw(10) << nrldata->points(i,0) << std::setw(10) << nrldata->points(i,1) << std::setw(10) << nrldata->points(i,2) << "\n";
        }
        for (unsigned int i=0; i<nrldata->nloads; ++i){
            std::cout << std::setw(10) << nrldata->boundaryconditions(i) << "\n";
        }
        for (unsigned int i=0; i<nrldata->npoints; ++i){
            std::cout << nrldata->exx(i+4*nrldata->npoints,5) << "\n";
        }
    }*/

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
}
