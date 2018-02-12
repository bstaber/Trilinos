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
    Epetra_SerialDenseMatrix points,exx,eyy,exy,energy;
    Epetra_SerialDenseVector boundaryconditions, angles;

    std::vector<int>    local_id_faces;
    std::vector<double> local_xi;
    std::vector<double> local_eta;

    distributenrldata(mesh & Mesh, std::string & path);
    ~distributenrldata();
    void retrieve_data(mesh & Mesh, std::string & path);
    double inverse_isoparametric_mapping(double & testx, double & testy, Epetra_SerialDenseVector & x, Epetra_SerialDenseVector & y, double & xi, double & eta);

};
#endif
