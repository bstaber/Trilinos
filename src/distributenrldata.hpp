#ifndef DISTRIBUTENRLDATA_HPP
#define DISTRIBUTENRLDATA_HPP

#include <stdio.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "readnrldata.hpp"
#include "meshpp.hpp"

class distributenrldata
{
public:
    
    int npoints,nloads;
    Epetra_SerialDenseMatrix points,exx,eyy,exy;
    Epetra_SerialDenseVector boundaryconditions, energy;
    
    std::vector<int> local_cells;
    std::vector<double> local_xi;
    std::vector<double> local_eta;

    distributenrldata(mesh & Mesh);
    ~distributenrldata();
    void retrieve_data(mesh & Mesh);
    double inverse_isoparametric_mapping(double & testx, double & testy, Epetra_SerialDenseVector & x, Epetra_SerialDenseVector & y, double & xi, double & eta);
    
};
#endif
