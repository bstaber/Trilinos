/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef NEWTONRAPHSON_HPP
#define NEWTONRAPHSON_HPP

#include "nonLinearFiniteElementProblem.hpp"
#include "Teuchos_RCP.hpp"
#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

class newtonRaphson
{
public:

    Epetra_Comm   * Comm;
    Epetra_Vector * x;

    nonLinearFiniteElementProblem * interface;
    Teuchos::ParameterList * Krylov;

    int MyPID, iter_max, iter_min, nb_bis_max;
    double time, delta, norm_inf_tol, norm_inf_max, eps, success, failure, umax, bc_disp, pressure_load, tol;

    newtonRaphson(nonLinearFiniteElementProblem & Interface, Teuchos::ParameterList & Parameters);
    ~newtonRaphson();

    int Solve_with_Aztec(bool print);
    //int Solve_with_Stratimikos(Teuchos::RCP<Teuchos::ParameterList> solverBuilderSL);
    int print_newton_solution(std::string fileName);
    void Initialization();
    void setInitialization(Epetra_Vector & init);
    void setParameters(Teuchos::ParameterList & Parameters);
};
#endif
