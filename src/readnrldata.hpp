/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef READNRLDATA_HPP
#define READNRLDATA_HPP

#include <stdio.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

class readnrldata
{
public:

    int npoints, nloads;
    Epetra_SerialDenseMatrix points, exx, eyy, exy, energy;
    Epetra_SerialDenseVector boundaryconditions, angles;

    readnrldata(bool load, std::string & path);
    ~readnrldata();

    void import_boundaryconditions(std::string & path);
    void import_exp_points(std::string & path);
    void import_expenergy(std::string & path);
    void import_exp_def(std::string & path);

};
#endif
