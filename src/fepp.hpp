/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef FEPP_HPP
#define FEPP_HPP

#include "meshpp.hpp"

namespace tri3{
    void shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta);
    void d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta);
    void dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & DX, double & jac, Epetra_SerialDenseMatrix & JacobianMatrix);
}

namespace quad4{
    void shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta);
    void d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta);
    void dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & DX, double & jac, Epetra_SerialDenseMatrix & JacobianMatrix);
}

namespace tri6{
    void shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta);
    void d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta);
    void dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & DX, double & jac, Epetra_SerialDenseMatrix & JacobianMatrix);
}

namespace tetra4{
    void shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta, double & zeta);
    void d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta, double & zeta);
}

namespace tetra10{
    void shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta, double & zeta);
    void d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta, double & zeta);
}

namespace hexa8{
    void shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta, double & zeta);
    void d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta, double & zeta);
}

void dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix JacobianMatrix, double & jac, Epetra_SerialDenseMatrix & DX);

void jacobian_matrix(Epetra_SerialDenseMatrix & X, Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & JacobianMatrix);
void jacobian_det(Epetra_SerialDenseMatrix & JacobianMatrix, double & jac);
void jacobian_det_faces(Epetra_SerialDenseMatrix & JacobianMatrix, double & jac);

void gauss_points_tri1(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta);
void gauss_points_tri3(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta);
void gauss_points_tri4(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta);

void gauss_points_quad1(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta);
void gauss_points_quad4(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta);
void gauss_points_quad9(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta);

void gauss_points_hexa4(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta);
void gauss_points_hexa8(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta);
void gauss_points_hexa27(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta);

void gauss_points_tetra1(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta);
void gauss_points_tetra4(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta);
void gauss_points_tetra5(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta);
void gauss_points_tetra11(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta);

int pnpoly(int & nvert, Epetra_SerialDenseVector & vertx, Epetra_SerialDenseVector & verty, double & testx, double & testy);


#endif
