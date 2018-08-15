#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "NRL_ModelF.hpp"
#include "newtonRaphson.hpp"
#include "objectiveFunction.hpp"
#include <random>

int main(int argc, char *argv[]){
    
    std::string    xmlInFileName = "";
    std::string    xmlExtraFileName = "";
    
    Teuchos::CommandLineProcessor  clp(false);
    clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
    clp.setOption("extra-file",&xmlExtraFileName,"Extra XML file to read into a parameter list");
    clp.setDocString("TO DO.");
    
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parse_return = clp.parse(argc,argv);
    if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;
        return parse_return;
    }
	
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    
    Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(new Teuchos::ParameterList);
    if(xmlInFileName.length()) {
        Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(*paramList));
    }
    if(xmlExtraFileName.length()) {
        Teuchos::updateParametersFromXmlFile(xmlExtraFileName, inoutArg(*paramList));
    }
    
    if (Comm.MyPID()==0){
        paramList->print(std::cout,2,true,true);
    }
    
    Teuchos::RCP<objectiveFunction<double>> obj = Teuchos::rcp(new objectiveFunction<double>(Comm,*paramList));
    
    if (Comm.MyPID()==0){
        int width = 20;
        int eglob = 4000;
        int nvert = obj->interface->Mesh->face_type;
        Epetra_SerialDenseVector x(nvert), y(nvert);
    
        double testx = 0.0;
        double testy = 0.0;
        double xmin, xmax;
        double ymin, ymax;
    
        int node;
        for (unsigned int inode=0; inode<nvert; ++inode){
            node = obj->interface->Mesh->faces_nodes[nvert*eglob+inode];
            x(inode) = obj->interface->Mesh->nodes_coord[3*node+0];
            y(inode) = obj->interface->Mesh->nodes_coord[3*node+1];
            testx += x(inode);
            testy += y(inode);
            if (inode==0){
                xmin = x(inode); xmax = x(inode);
                ymin = y(inode); ymax = y(inode);
            }
            else{
                if (x(inode)<xmin){
                    xmin = x(inode);
                }
                if (x(inode)>xmax){
                    xmax = x(inode);
                }
                if (y(inode)<ymin){
                    ymin = y(inode);
                }
                if (y(inode)>ymax){
                    ymax = y(inode);
                }
            }
        }
        
        std::mt19937 G(time(NULL));
        std::uniform_real_distribution<double>  u(xmin,xmax);
        std::uniform_real_distribution<double>  v(ymin,ymax);
        
        std::cout << std::setw(width) << "cell" << std::setw(width) << "testx" << std::setw(width) << "testy" << std::setw(width) << "pnpoly" << std::setw(width) << "xi" << std::setw(width) << "eta" << std::setw(width) << "residual\n";
        int nnmc = 500;
        for (unsigned int nmc=0; nmc<nnmc; ++nmc){
            testx = u(G);
            testy = v(G);
        
            int result = pnpoly(nvert,x,y,testx,testy);
            /*if (!result){
                std::cout << "**ERR: I'm not in the cell!\n";
             }*/
    
            double xi = 0.0;
            double eta = 0.0;
            double residual = 0.0;
            if (result){
                residual = obj->inverse_isoparametric_mapping(testx, testy, x, y, xi, eta);
            }
    
            std::cout << std::setw(width) << eglob << std::setw(width) << testx << std::setw(width) << testy << std::setw(width) << result << std::setw(width) << xi << std::setw(width) << eta << std::setw(width) << residual << "\n";
        }
    }
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;
    
}
