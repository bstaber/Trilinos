#ifndef NEWTON_RAPHSONPP_HPP
#define NEWTON_RAPHSONPP_HPP

#include "hyperelasticity_setup_pp.hpp"
#include "Finite_Element_Problem.hpp"
#include "Teuchos_RCP.hpp"
#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosEpetraAdapter.hpp"
#include <BelosSolverFactory.hpp>
#include "BelosBlockGmresSolMgr.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

class Newton_Raphson
{
public:
    Newton_Raphson(Finite_Element_Problem & Interface, Teuchos::ParameterList & Parameters);
    ~Newton_Raphson();
    
    int Solve_with_Aztec();
    int Solve_with_Stratimikos(Teuchos::RCP<Teuchos::ParameterList> solverBuilderSL);
    int print_newton_solution(std::string fileName);
    void Initialization();
    void setInitialization(Epetra_Vector & init);
    void setParameters(Teuchos::ParameterList & Parameters);
        
    int iter_max, iter_min;
    int nb_bis_max;
    int MyPID;
    double time;
    double delta;
    double norm_inf_tol;
    double norm_inf_max;
    double eps;
    double success;
    double failure;
    double umax;
    double bc_disp;
    double pressure_load;
    
    Epetra_Comm * Comm;
    Epetra_Vector * x;
    
    Finite_Element_Problem * interface;
    Teuchos::ParameterList * Krylov;
};
#endif
