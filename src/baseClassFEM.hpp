/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef BASECLASSFEM_HPP
#define BASECLASSFEM_HPP
#include "meshpp.hpp"
#include "error_messager.hpp"
class baseClassFEM
{
public:
    baseClassFEM();
    ~baseClassFEM();

    mesh * Mesh;
    Epetra_Comm * Comm;

    Epetra_Map * OverlapMap;
    Epetra_Map * StandardMap;
    Epetra_Import * ImportToOverlapMap;
    Epetra_FECrsGraph * FEGraph;
};
#endif
