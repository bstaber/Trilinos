#ifndef NRL_MODELF_HPP
#define NRL_MODELF_HPP

#include "tensor_calculus.hpp"
#include "hyperelasticity_setup_pp.hpp"

class NRL_ModelF : public hyperelasticity_setup
{
public:
    
    NRL_ModelF(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){
        
    }
    
    ~NRL_ModelF(){
    }
    
    void set_parameters(Epetra_SerialDenseVector & x){
        
    }
    
    void set_plyagl(double & Plyagl){
        
    }
    
    void get_matrix_and_rhs(Epetra_Vector & x, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
        
    }
    
    void setup_dirichlet_conditions(){
        
    }
    
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
        
    }
    
    void get_material_parameters(unsigned int & e_lid, unsigned int & gp){
        
    }
    
    void get_constitutive_tensors(Epetra_SerialDenseMatrix & deformation_gradient, Epetra_SerialDenseVector &
                                  
    }
    
    void get_constitutive_tensors_static_condensation(Epetra_SerialDenseMatrix & deformation_gradient, double & det,
    }
    
    void get_internal_pressure(double & theta, double & pressure, double & dpressure){
        std::cerr << "**Err: Not using static condensation method!\n";
    }
    
    void get_material_parameters_for_recover(unsigned int & e_lid, double & xi, double & eta, double & zeta){
    }
                                                      
    void get_stress_for_recover(Epetra_SerialDenseMatrix & deformation_gradient, double & det, Epetra_SerialDenseMatrix & piola_stress){
    }
    
    double pr;
    double m1;
    double m2;
    double m;
    double beta3;
    double beta4;
    double beta5;
    
    double pmbeta4;
    double pmbeta5;
    
    double plyagl;
    double cos_plyagl;
    double sin_plyagl;
    int n_ply;
    std::vector<int> phase;
    
};

#endif
