#ifndef LINEAR_FINITE_ELEMENT_PROBLEM_HPP
#define LINEAR_FINITE_ELEMENT_PROBLEM_HPP
#include "BaseClassFEM.hpp"
class Linear_Finite_Element_Problem : public BaseClassFEM
{
public:
    Linear_Finite_Element_Problem(){
    };
    ~Linear_Finite_Element_Problem(){
    };
    virtual void setup_dirichlet_conditions(){
    };
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
    };
};
#endif
