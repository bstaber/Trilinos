#include "readnrldata.hpp"
#include <math.h>

readnrldata::readnrldata(){
    import_boundaryconditions();
    import_exp_points();
    import_exp_def();
}

readnrldata::~readnrldata(){
}

void readnrldata::import_boundaryconditions(){
    
    std::string filename = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/dirichletbcs.txt";
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
}

void readnrldata::import_exp_points(){
    double x,y,z;
    std::string filename = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/xyz.txt";
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

void readnrldata::import_exp_def(){
    
    std::string path = "/Users/brian/Documents/GitHub/Trilinos_results/nrl/data/";
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


