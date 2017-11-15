#ifndef BASECLASSFEM_HPP
#define BASECLASSFEM_HPP

#include "meshpp.hpp"

class BaseClassFEM
{
    
public:
    BaseClassFEM(){
    };
    ~BaseClassFEM(){
    };
    
    mesh * Mesh;
    Epetra_Comm * Comm;
    
    Epetra_Map * OverlapMap;
    Epetra_Map * StandardMap;
    Epetra_Import * ImportToOverlapMap;
    Epetra_FECrsGraph * FEGraph;
};

#endif
