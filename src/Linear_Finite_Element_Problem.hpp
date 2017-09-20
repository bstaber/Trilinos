#ifndef LINEAR_FINITE_ELEMENT_PROBLEM_HPP
#define LINEAR_FINITE_ELEMENT_PROBLEM_HPP

#include "meshpp.hpp"

class Linear_Finite_Element_Problem
{
    
public:
    Linear_Finite_Element_Problem(){
    };
    ~Linear_Finite_Element_Problem(){
    };
    
    Epetra_SerialDenseVector dead_pressure;
        
    mesh * Mesh;
    Epetra_Comm * Comm;
    
    Epetra_Map * OverlapMap;
    Epetra_Map * StandardMap;
    Epetra_Import * ImportToOverlapMap;
    Epetra_FECrsGraph * FEGraph;
};

#endif
