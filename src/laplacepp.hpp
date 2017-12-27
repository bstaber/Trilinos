#ifndef LAPLACEPP_HPP
#define LAPLACEPP_HPP

#include "Linear_Finite_Element_Problem.hpp"

class laplace : public Linear_Finite_Element_Problem
{
public:
    laplace();
    laplace(mesh & Mesh);
    laplace(mesh & Mesh, Epetra_SerialDenseMatrix & bc_matrix);
    laplace(mesh & Mesh, Teuchos::ParameterList & Parameters);
    laplace(Epetra_Comm & comm, Teuchos::ParameterList & Parameters);
    ~laplace();
    void create_FECrsGraph();
    void assembling(Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs);
    void assembling_OAZ(Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs, int * bc_indx, double * bc_val);
    void solve_aztec(Teuchos::ParameterList & Parameters,  Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs, int * bc_indx, double * bc_val);
    void solve_amesos(Teuchos::ParameterList & Parameters, Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs, int * bc_indx, double * bc_val);
    void compute_local_directions(Epetra_Vector & laplace_one, Epetra_Vector & laplace_two);
    void compute_center_local_directions(Epetra_Vector & laplace_one, Epetra_Vector & laplace_two);
    int print_solution(Epetra_Vector & lhs, std::string fileName);

    Epetra_SerialDenseMatrix laplace_direction_one, laplace_direction_one_center;
    Epetra_SerialDenseMatrix laplace_direction_two, laplace_direction_two_center;
    Epetra_SerialDenseMatrix laplace_direction_two_cross_one, laplace_direction_two_cross_one_center;
};
#endif
