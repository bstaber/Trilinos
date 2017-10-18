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
    
    int npoints,nloads;
    Epetra_SerialDenseMatrix points,exx,eyy,exy;
    Epetra_SerialDenseVector boundaryconditions,energy;
    
    readnrldata();
    ~readnrldata();
    
    void import_boundaryconditions();
    void import_exp_points();
    void import_expenergy();
    void import_exp_def();
    
};
#endif
