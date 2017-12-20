#include "readnrldata.hpp"
#include <math.h>

readnrldata::readnrldata(bool load, std::string & path){
    if(load){
        import_boundaryconditions(path);
        import_exp_points(path);
        import_expenergy(path);
        import_exp_def(path);
        angles.Resize(8);
        angles(0) = 30.0;
        angles(1) = 60.0;
        angles(2) = 60.0;
        angles(3) = 30.0;
        angles(4) = 15.0;
        angles(5) = 15.0;
        angles(6) = 75.0;
        angles(7) = 75.0;
    }
}

readnrldata::~readnrldata(){
}

void readnrldata::import_boundaryconditions(std::string & path){
    //std::string filename = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/dirichletbcs.txt";
    //std::string filename = "/home/s/staber/Trilinos_results/nrl/data/dirichletbcs.txt";
    std::string filename = path + "dirichletbcs.txt";
    std::ifstream file;
    double gbc;
    file.open(filename);
    if (file.is_open()){
        file >> nloads;
        boundaryconditions.Resize(nloads);
        for (unsigned int i=0; i<nloads; ++i){
            file >> gbc;
            boundaryconditions(i) = gbc;
        }
        file.close();
    }
    else{
        std::cout << "Couldn't open the file containing the boundary conditions.\n";
    }
    
/*    filename = path + "linspace.txt";
    file.open(filename);
    if (file.is_open()){
        linspace.Resize(nloads);
        for (unsigned int i=0; i<nloads; ++i){
            file >> gbc;
            linspace(i) = gbc;
        }
    }
*/
}

void readnrldata::import_expenergy(std::string & path){
    //std::string path = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/expenergy.txt";
    //std::string path = "/home/s/staber/Trilinos_results/nrl/data/expenergy.txt";
    std::string filename = path + "expenergy.txt";
    std::ifstream file;
    double gen;
    file.open(filename);
    if (file.is_open()){
        energy.Reshape(nloads,8);
        for (unsigned int i=0; i<nloads; ++i){
            for (unsigned int j=0; j<8; ++j){
                file >> gen;
                energy(i,j) = gen;
            }
        }
        file.close();
    }
    else{
        std::cout << "Couldn't open the file containing the experimental energies.\n";
    }
}

void readnrldata::import_exp_points(std::string & path){
    double x,y,z;
    //std::string filename = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/xyz.txt";
    //std::string filename = "/home/s/staber/Trilinos_results/nrl/data/xyz.txt";
    std::string filename = path + "xyz.txt";
    std::ifstream file;
    file.open(filename);
    if (file.is_open()){
        file >> npoints;
        points.Reshape(npoints,3);
        for (unsigned int i=0; i<npoints; ++i){
            file >> x;
            file >> y;
            file >> z;
            points(i,0) = x;
            points(i,1) = y;
            points(i,2) = z;
        }
        file.close();
    }
    else{
        std::cout << "Couldn't open the file containing the locations of the experimental points.\n";
    }
}

void readnrldata::import_exp_def(std::string & path){
    //std::string path = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/";
    //std::string path = "/home/s/staber/Trilinos_results/nrl/data/";
    double gexx,geyy,gexy;
    for (unsigned int id=0; id<8; ++id){
        std::string path_exx = path + "exx_id" + std::to_string(id+1) + ".txt";
        std::string path_eyy = path + "eyy_id" + std::to_string(id+1) + ".txt";
        std::string path_exy = path + "exy_id" + std::to_string(id+1) + ".txt";
        std::ifstream file_exx,file_eyy,file_exy;
        file_exx.open(path_exx); file_eyy.open(path_eyy); file_exy.open(path_exy);
        exx.Reshape(nloads*npoints,8); eyy.Reshape(nloads*npoints,8); exy.Reshape(nloads*npoints,8);
        if (file_exx.is_open() && file_eyy.is_open() && file_exy.is_open()){
            for (unsigned int t=0; t<nloads; ++t){
                for (unsigned int j=0; j<npoints; ++j){
                    file_exx >> gexx; file_eyy >> geyy; file_exy >> gexy;
                    exx(j+t*npoints,id) = gexx; eyy(j+t*npoints,id) = geyy; exy(j+t*npoints,id) = gexy;
                }
            }
            
        }
        else{
            std::cout << "Couldn't open one or some of the eij files.\n";
        }
    }
}


