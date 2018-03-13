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

    //std::string mesh_file = "/Users/brian/Documents/GitHub/Trilinos/cee530/mesh/manufactured1.msh";
    std::string mesh_file = "/Users/brian/Documents/GitHub/Trilinos/arteries/mesh/media_flatboundaries.msh";
    //std::string mesh_file = "/Users/brian/Documents/GitHub/Trilinos/nrl/mesh/composite_hexa_32.msh";
    mesh Mesh(Comm, mesh_file, 1.0);

    Comm.Barrier();
    if (Comm.MyPID()==0){
        std::cout << "EPART:\n";
        for (unsigned int i=0; i<Mesh.n_cells; ++i){
            std::cout << Mesh.epart[i] << "\n";
        }
        std::cout << "FACES:\n";
        for (unsigned int i=0; i<Mesh.n_local_faces; ++i){
            std::cout << Mesh.local_faces[i] << "\n";
        }
    }

    /*std::vector<int> phase;
    for (unsigned int e=0; e<Mesh.n_cells/32; ++e){
        for (unsigned int j=0; j<32; ++j){
            phase.push_back(j);
        }
    }

    for (unsigned int e_gid=0; e_gid<32; ++e_gid){
        if (phase[e_gid] % 2){
            std::cout << e_gid << std::setw(15) << phase[e_gid] << std::setw(15) << "true\n";
        }
        else{
            std::cout << e_gid << std::setw(15) << phase[e_gid] << std::setw(15) << "false\n";
        }
    }*/

#ifdef HAVE_MPI
MPI_Finalize();
#endif
return 0;

}
