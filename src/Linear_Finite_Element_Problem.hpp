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
    
    //double pressure_load;
    Epetra_SerialDenseVector dead_pressure;
    
    virtual void setup_dirichlet_conditions() = 0;
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement) = 0;
    
    mesh * Mesh;
    Epetra_Comm * Comm;
    
    Epetra_Map * OverlapMap;
    Epetra_Map * StandardMap;
    Epetra_Import * ImportToOverlapMap;
    Epetra_FECrsGraph * FEGraph;
};

#endif
