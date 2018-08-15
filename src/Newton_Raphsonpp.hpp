#ifndef NEWTON_RAPHSONPP_HPP
#define NEWTON_RAPHSONPP_HPP

#include "nonLinearFiniteElementProblem.hpp"
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

    Epetra_Comm   * Comm;
    Epetra_Vector * x;

    nonLinearFiniteElementProblem * interface;
    Teuchos::ParameterList * Krylov;

    int MyPID, iter_max, iter_min, nb_bis_max;
    double time, delta, norm_inf_tol, norm_inf_max, eps, success, failure, umax, bc_disp, pressure_load, tol;

    Newton_Raphson(bonLinearFiniteElementProblem & Interface, Teuchos::ParameterList & Parameters);
    ~Newton_Raphson();

    int Solve_with_Aztec(bool print);
    int Solve_with_Stratimikos(Teuchos::RCP<Teuchos::ParameterList> solverBuilderSL);
    int print_newton_solution(std::string fileName);
    void Initialization();
    void setInitialization(Epetra_Vector & init);
    void setParameters(Teuchos::ParameterList & Parameters);
};
#endif
