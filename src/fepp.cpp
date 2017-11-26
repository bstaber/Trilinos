#include "fepp.hpp"

void tri3::shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta){
    N(0) = 1.0 - xi - eta;
    N(1) = xi;
    N(2) = eta;
}
void tri3::d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta){
    D(0,0) = -1.0; D(0,1) = -1.0;
    D(1,0) = 1.0;  D(1,1) = 0.0;
    D(2,0) = 0.0;  D(2,1) = 1.0;
}
void tri3::dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & DX, double & jac, Epetra_SerialDenseMatrix & JacobianMatrix){
    
    Epetra_SerialDenseMatrix InverseJacobianMatrix(2,2);
    InverseJacobianMatrix(0,0) = (1.0/jac)*JacobianMatrix(1,1);
    InverseJacobianMatrix(1,1) = (1.0/jac)*JacobianMatrix(0,0);
    InverseJacobianMatrix(0,1) = -(1.0/jac)*JacobianMatrix(0,1);
    InverseJacobianMatrix(1,0) = -(1.0/jac)*JacobianMatrix(1,0);
    
    DX.Multiply('N','N',1.0,D,InverseJacobianMatrix,0.0);
}

void quad4::shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta){
    N(0) = (1.0/4.0)*(1.0-xi)*(1.0-eta);
    N(1) = (1.0/4.0)*(1.0+xi)*(1.0-eta);
    N(2) = (1.0/4.0)*(1.0+xi)*(1.0+eta);
    N(3) = (1.0/4.0)*(1.0-xi)*(1.0+eta);
}
void quad4::d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta){
    D(0,0) = -(1.0/4.0)*(1.0-eta); D(0,1) = -(1.0/4.0)*(1.0-xi);
    D(1,0) = (1.0/4.0)*(1.0-eta);  D(1,1) = -(1.0/4.0)*(1.0+xi);
    D(2,0) = (1.0/4.0)*(1.0+eta);  D(2,1) = (1.0/4.0)*(1.0+xi);
    D(3,0) = -(1.0/4.0)*(1.0+eta); D(3,1) = (1.0/4.0)*(1.0-xi);
}
void quad4::dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & DX, double & jac, Epetra_SerialDenseMatrix & JacobianMatrix){
    
    Epetra_SerialDenseMatrix InverseJacobianMatrix(2,2);
    InverseJacobianMatrix(0,0) = (1.0/jac)*JacobianMatrix(1,1);
    InverseJacobianMatrix(1,1) = (1.0/jac)*JacobianMatrix(0,0);
    InverseJacobianMatrix(0,1) = -(1.0/jac)*JacobianMatrix(0,1);
    InverseJacobianMatrix(1,0) = -(1.0/jac)*JacobianMatrix(1,0);
    
    DX.Multiply('N','N',1.0,D,InverseJacobianMatrix,0.0);
}

void tri6::shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta){
    double lambda = 1.0 - xi - eta;
    N(0) = lambda*(2.0*lambda-1.0);
    N(1) = xi*(2.0*xi-1.0);
    N(2) = eta*(2.0*eta-1.0);
    N(3) = 4.0*lambda*xi;
    N(4) = 4.0*eta*xi;
    N(5) = 4.0*lambda*eta;
}

void tri6::d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta){
    D(0,0) = 4.0*xi + 4.0*eta - 3.0;
    D(1,0) = 4.0*xi - 1.0;
    D(2,0) = 0.0;
    D(3,0) = 4.0 - 8.0*xi - 4.0*eta;
    D(4,0) = 4.0*eta;
    D(5,0) = -4.0*eta;
    
    D(0,1) = 4.0*eta + 4.0*xi - 3.0;
    D(1,1) = 0.0;
    D(2,1) = 4.0*eta - 1.0;
    D(3,1) = -4.0*xi;
    D(4,1) = 4.0*xi;
    D(5,1) = 4.0 - 4.0*xi - 8.0*eta;
}

void tri6::dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & DX, double & jac, Epetra_SerialDenseMatrix & JacobianMatrix){
    
    Epetra_SerialDenseMatrix InverseJacobianMatrix(2,2);
    InverseJacobianMatrix(0,0) = (1.0/jac)*JacobianMatrix(1,1);
    InverseJacobianMatrix(1,1) = (1.0/jac)*JacobianMatrix(0,0);
    InverseJacobianMatrix(0,1) = -(1.0/jac)*JacobianMatrix(0,1);
    InverseJacobianMatrix(1,0) = -(1.0/jac)*JacobianMatrix(1,0);
    
    DX.Multiply('N','N',1.0,D,InverseJacobianMatrix,0.0);
}

void hexa8::shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta, double & zeta){
    N(0) = (1.0/8.0)*(1.0-xi)*(1.0-eta)*(1.0-zeta);
    N(1) = (1.0/8.0)*(1.0+xi)*(1.0-eta)*(1.0-zeta);
    N(2) = (1.0/8.0)*(1.0+xi)*(1.0+eta)*(1.0-zeta);
    N(3) = (1.0/8.0)*(1.0-xi)*(1.0+eta)*(1.0-zeta);
    N(4) = (1.0/8.0)*(1.0-xi)*(1.0-eta)*(1.0+zeta);
    N(5) = (1.0/8.0)*(1.0+xi)*(1.0-eta)*(1.0+zeta);
    N(6) = (1.0/8.0)*(1.0+xi)*(1.0+eta)*(1.0+zeta);
    N(7) = (1.0/8.0)*(1.0-xi)*(1.0+eta)*(1.0+zeta);
}

void hexa8::d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta, double & zeta){
    D(0,0) = -(1.0/8.0)*(1.0-eta)*(1.0-zeta); D(0,1) = -(1.0/8.0)*(1.0-xi)*(1.0-zeta); D(0,2) = -(1.0/8.0)*(1.0-xi)*(1.0-eta);
    D(1,0) = (1.0/8.0)*(1.0-eta)*(1.0-zeta);  D(1,1) = -(1.0/8.0)*(1.0+xi)*(1.0-zeta); D(1,2) = -(1.0/8.0)*(1.0+xi)*(1.0-eta);
    D(2,0) = (1.0/8.0)*(1.0+eta)*(1.0-zeta);  D(2,1) = (1.0/8.0)*(1.0+xi)*(1.0-zeta);  D(2,2) = -(1.0/8.0)*(1.0+xi)*(1.0+eta);
    D(3,0) = -(1.0/8.0)*(1.0+eta)*(1.0-zeta); D(3,1) = (1.0/8.0)*(1.0-xi)*(1.0-zeta);  D(3,2) = -(1.0/8.0)*(1.0-xi)*(1.0+eta);
    D(4,0) = -(1.0/8.0)*(1.0-eta)*(1.0+zeta); D(4,1) = -(1.0/8.0)*(1.0-xi)*(1.0+zeta); D(4,2) = (1.0/8.0)*(1.0-xi)*(1.0-eta);
    D(5,0) = (1.0/8.0)*(1.0-eta)*(1.0+zeta);  D(5,1) = -(1.0/8.0)*(1.0+xi)*(1.0+zeta); D(5,2) = (1.0/8.0)*(1.0+xi)*(1.0-eta);
    D(6,0) = (1.0/8.0)*(1.0+eta)*(1.0+zeta);  D(6,1) = (1.0/8.0)*(1.0+xi)*(1.0+zeta);  D(6,2) = (1.0/8.0)*(1.0+xi)*(1.0+eta);
    D(7,0) = -(1.0/8.0)*(1.0+eta)*(1.0+zeta); D(7,1) = (1.0/8.0)*(1.0-xi)*(1.0+zeta);  D(7,2) = (1.0/8.0)*(1.0-xi)*(1.0+eta);
}

void tetra4::shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta, double & zeta){
    double lambda = 1 - xi - eta - zeta;
    N(0) = lambda;
    N(1) = xi;
    N(2) = eta;
    N(3) = zeta;
}

void tetra4::d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta, double & zeta){
    D(0,0) = -1.0; D(0,1) = -1.0; D(0,2) = -1.0;
    D(1,0) = 1.0; D(1,1) = 0.0; D(1,2) = 0.0;
    D(2,0) = 0.0; D(2,1) = 1.0; D(2,2) = 0.0;
    D(3,0) = 0.0; D(3,1) = 0.0; D(3,2) = 1.0;
}

void tetra10::shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta, double & zeta){
    double lambda = 1 - xi - eta - zeta;
    N(0) = lambda*(2.0*lambda-1.0);
    N(1) = xi*(2.0*xi-1.0);
    N(2) = eta*(2.0*eta-1.0);
    N(3) = zeta*(2.0*zeta-1.0);
    N(4) = 4*xi*lambda;
    N(5) = 4*xi*eta;
    N(6) = 4*eta*lambda;
    N(7) = 4*zeta*lambda;
    N(8) = 4*eta*zeta;
    N(9) = 4*xi*zeta;
}

void tetra10::d_shape_functions(Epetra_SerialDenseMatrix & D, double & xi, double & eta, double & zeta){
    D(0,0) = 4.0*eta + 4.0*xi + 4.0*zeta - 3.0; D(0,1) = 4.0*eta + 4.0*xi + 4.0*zeta - 3.0;  D(0,2) = 4.0*eta + 4.0*xi + 4.0*zeta - 3.0;
    D(1,0) = 4.0*xi - 1.0; D(1,1) = 0.0; D(1,2) = 0.0;
    D(2,0) = 0.0; D(2,1) = 4.0*eta - 1.0; D(2,2) = 0.0;
    D(3,0) = 0.0; D(3,1) = 0.0; D(3,2) = 4.0*zeta - 1.0;
    D(4,0) = 4.0 - 8.0*xi - 4.0*zeta - 4.0*eta; D(4,1) = -4.0*xi; D(4,2) = -4.0*xi;
    D(5,0) = 4.0*eta; D(5,1) = 4.0*xi; D(5,2) = 0.0;
    D(6,0) = -4.0*eta; D(6,1) = 4.0 - 4.0*xi - 4.0*zeta - 8.0*eta; D(6,2) = -4.0*eta;
    D(7,0) = -4.0*zeta; D(7,1) = -4.0*zeta; D(7,2) = 4.0 - 4.0*xi - 8.0*zeta - 4.0*eta;
    D(8,0) = 0.0; D(8,1) = 4.0*zeta; D(8,2) = 4.0*eta;
    D(9,0) = 4.0*zeta; D(9,1) = 0.0; D(9,2) = 4.0*xi;
}

void dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix JacobianMatrix, double & jac, Epetra_SerialDenseMatrix & DX)
{
    Epetra_SerialDenseMatrix InverseJacobianMatrix(3,3);
    
    InverseJacobianMatrix(0,0) = (1.0/jac)*(JacobianMatrix(1,1)*JacobianMatrix(2,2)-JacobianMatrix(1,2)*JacobianMatrix(2,1));
    InverseJacobianMatrix(0,1) = (1.0/jac)*(JacobianMatrix(0,2)*JacobianMatrix(2,1)-JacobianMatrix(0,1)*JacobianMatrix(2,2));
    InverseJacobianMatrix(0,2) = (1.0/jac)*(JacobianMatrix(0,1)*JacobianMatrix(1,2)-JacobianMatrix(0,2)*JacobianMatrix(1,1));
    
    InverseJacobianMatrix(1,0) = (1.0/jac)*(JacobianMatrix(1,2)*JacobianMatrix(2,0)-JacobianMatrix(1,0)*JacobianMatrix(2,2));
    InverseJacobianMatrix(1,1) = (1.0/jac)*(JacobianMatrix(0,0)*JacobianMatrix(2,2)-JacobianMatrix(0,2)*JacobianMatrix(2,0));
    InverseJacobianMatrix(1,2) = (1.0/jac)*(JacobianMatrix(0,2)*JacobianMatrix(1,0)-JacobianMatrix(0,0)*JacobianMatrix(1,2));
    
    InverseJacobianMatrix(2,0) = (1.0/jac)*(JacobianMatrix(1,0)*JacobianMatrix(2,1)-JacobianMatrix(1,1)*JacobianMatrix(2,0));
    InverseJacobianMatrix(2,1) = (1.0/jac)*(JacobianMatrix(0,1)*JacobianMatrix(2,0)-JacobianMatrix(0,0)*JacobianMatrix(2,1));
    InverseJacobianMatrix(2,2) = (1.0/jac)*(JacobianMatrix(0,0)*JacobianMatrix(1,1)-JacobianMatrix(0,1)*JacobianMatrix(1,0));
    
    DX.Multiply('N','N',1.0,D,InverseJacobianMatrix,0.0);
}

void jacobian_matrix(Epetra_SerialDenseMatrix & X, Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & JacobianMatrix){
JacobianMatrix.Multiply('N','N',1.0,X,D,0.0);
}

void jacobian_det(Epetra_SerialDenseMatrix & JacobianMatrix, double & jac){
    jac = fabs(JacobianMatrix(0,0)*JacobianMatrix(1,1)*JacobianMatrix(2,2)-JacobianMatrix(0,0)*JacobianMatrix(1,2)*JacobianMatrix(2,1)-JacobianMatrix(0,1)*JacobianMatrix(1,0)*JacobianMatrix(2,2)+JacobianMatrix(0,1)*JacobianMatrix(1,2)*JacobianMatrix(2,0)+JacobianMatrix(0,2)*JacobianMatrix(1,0)*JacobianMatrix(2,1)-JacobianMatrix(0,2)*JacobianMatrix(1,1)*JacobianMatrix(2,0) );
}

void jacobian_det_faces(Epetra_SerialDenseMatrix & JacobianMatrix, double & jac){
    jac = fabs(JacobianMatrix(0,0)*JacobianMatrix(1,1) - JacobianMatrix(0,1)*JacobianMatrix(0,1));
}

void gauss_points_tri1(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta){
    weight.Resize(1); xi.Resize(1); eta.Resize(1);
    weight(0) = 1.0/2.0;
    xi(0) = 1.0/3.0;
    eta(0) = 1.0/3.0;
}
void gauss_points_tri3(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta){
    weight.Resize(3); xi.Resize(3); eta.Resize(3);
    weight(0) = 1.0/6.0; weight(1) = 1.0/6.0; weight(2) = 1.0/6.0;
    xi(0) = 1.0/6.0; xi(1) = 2.0/3.0; xi(2) = 1.0/6.0;
    eta(0) = 1.0/6.0; eta(1) = 1.0/6.0; eta(2) = 2.0/3.0;
}
void gauss_points_tri4(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta){
    weight.Resize(4); xi.Resize(4); eta.Resize(4);
    weight(0) = -27.0/96.0; weight(1) = 25.0/96.0; weight(2) = 25.0/96.0; weight(3) = 25.0/96.0;
    xi(0) = 1.0/3.0; xi(1) = 1.0/5.0; xi(2) = 3.0/5.0; xi(3) = 1.0/5.0;
    eta(0) = 1.0/3.0; eta(1) = 1.0/5.0; eta(2) = 1.0/5.0; eta(3) = 3.0/5.0;
}

void gauss_points_quad1(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta){
    weight.Resize(1); xi.Resize(1); eta.Resize(1);
    weight(0) = 4.0;
    xi(0) = 0.0;
    eta(0) = 0.0;
}
void gauss_points_quad4(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta){
    weight.Resize(4); xi.Resize(4); eta.Resize(4);
    double a = 1.0/std::sqrt(3.0);
    weight(0) = 1.0; weight(1) = 1.0; weight(2) = 1.0; weight(3) = 1.0;
    xi(0) = -a; xi(1) = a; xi(2) = -a; xi(3) = a;
    eta(0) = -a; eta(1) = -a; eta(2) = a; eta(3) = a;
}
void gauss_points_quad9(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta){
    weight.Resize(9); xi.Resize(9); eta.Resize(9);
    double a = std::sqrt(3.0/5.0);
    weight(0) = 25.0/81.0; weight(1) = 40.0/81.0; weight(2) = 25.0/81.0; weight(3) = 40.0/81.0; weight(4) = 64.0/81.0; weight(5) = 40.0/81.0; weight(6) = 25.0/81.0; weight(7) = 40.0/81.0; weight(8) = 25.0/81.0;
    xi(0) = -a; xi(1) = 0.0; xi(2) = a; xi(3) = -a; xi(4) = 0.0; xi(5) = a; xi(6) = -a; xi(7) = 0.0; xi(8) = a;
    eta(0) = -a; eta(1) = -a; eta(2) = -a; eta(3) = 0.0; eta(4) = 0.0; eta(5) = 0.0; eta(6) = a; eta(7) = a; eta(8) = a;
}

void gauss_points_hexa4(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta){
    weight.Resize(4); xi.Resize(4); eta.Resize(4); zeta.Resize(4);
    double alpha = std::sqrt(2.0/3.0);
    double beta = std::sqrt(1.0/3.0);
    weight(0) = 2.0; weight(1) = 2.0; weight(2) = 2.0; weight(3) = 2.0;
    xi(0) = 0.0; eta(0) = alpha; zeta(0) = -beta;
    xi(1) = 0.0; eta(1) = -alpha; zeta(1) = -beta;
    xi(2) = alpha; eta(2) = 0.0; zeta(2) = beta;
    xi(3)= -alpha; eta(3) = 0.0; zeta(3) = beta;
}
void gauss_points_hexa6(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta){
    weight.Resize(6); xi.Resize(6); eta.Resize(6); zeta.Resize(6);
    //TODO
}
void gauss_points_hexa8(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta){
    weight.Resize(8); xi.Resize(8); eta.Resize(8); zeta.Resize(8);
    
    double a = 1.0/std::sqrt(3.0);
    xi(0)=-a;  eta(0)=-a;  zeta(0)=-a;   weight(0)=1.0;
    xi(1)=-a;  eta(1)=-a;  zeta(1)=a;    weight(1)=1.0;
    xi(2)=-a;  eta(2)=a;   zeta(2)=-a;   weight(2)=1.0;
    xi(3)=-a;  eta(3)=a;   zeta(3)=a;    weight(3)=1.0;
    xi(4)= a;  eta(4)=-a;  zeta(4)=-a;   weight(4)=1.0;
    xi(5)= a;  eta(5)=-a;  zeta(5)=a;    weight(5)=1.0;
    xi(6)= a;  eta(6)=a;   zeta(6)=-a;   weight(6)=1.0;
    xi(7)= a;  eta(7)=a;   zeta(7)=a;    weight(7)=1.0;
}

void gauss_points_hexa27(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta){
    weight.Resize(27); xi.Resize(27); eta.Resize(27); zeta.Resize(27);
    
    double a = std::sqrt(3.0/5.0); double c1 = 5.0/9.0; double c2 = 8.0/9.0;
    xi(0) =   -a;  eta(0) =  -a;  zeta(0) =  -a;  weight(0) = c1*c1*c1;
    xi(1) =   -a;  eta(1) =  -a;  zeta(1) = 0.0;  weight(1) = c1*c1*c2;
    xi(2) =   -a;  eta(2) =  -a;  zeta(2) =   a;  weight(2) = c1*c1*c1;
    xi(3) =   -a;  eta(3) = 0.0;  zeta(3) =  -a;  weight(3) = c1*c1*c2;
    xi(4) =   -a;  eta(4) = 0.0;  zeta(4) = 0.0;  weight(4) = c1*c2*c2;
    xi(5) =   -a;  eta(5) = 0.0;  zeta(5) =   a;  weight(5) = c1*c1*c2;
    xi(6) =   -a;  eta(6) =   a;  zeta(6) =  -a;  weight(6) = c1*c1*c1;
    xi(7) =   -a;  eta(7) =   a;  zeta(7) = 0.0;  weight(7) = c1*c1*c2;
    xi(8) =   -a;  eta(8) =   a;  zeta(8) =   a;  weight(8) = c1*c1*c1;
    xi(9) =  0.0;  eta(9) =  -a;  zeta(9) =  -a;  weight(9) = c1*c1*c2;
    xi(10) = 0.0; eta(10) =  -a; zeta(10) = 0.0; weight(10) = c1*c2*c2;
    xi(11) = 0.0; eta(11) =  -a; zeta(11) =   a; weight(11) = c1*c1*c2;
    xi(12) = 0.0; eta(12) = 0.0; zeta(12) =  -a; weight(12) = c1*c2*c2;
    xi(13) = 0.0; eta(13) = 0.0; zeta(13) = 0.0; weight(13) = c2*c2*c2;
    xi(14) = 0.0; eta(14) = 0.0; zeta(14) =   a; weight(14) = c1*c2*c2;
    xi(15) = 0.0; eta(15) =   a; zeta(15) =  -a; weight(15) = c1*c1*c2;
    xi(16) = 0.0; eta(16) =   a; zeta(16) = 0.0; weight(16) = c1*c2*c2;
    xi(17) = 0.0; eta(17) =   a; zeta(17) =   a; weight(17) = c1*c1*c2;
    xi(18) =   a; eta(18) =  -a; zeta(18) =  -a; weight(18) = c1*c1*c1;
    xi(19) =   a; eta(19) =  -a; zeta(19) = 0.0; weight(19) = c1*c1*c2;
    xi(20) =   a; eta(20) =  -a; zeta(20) =   a; weight(20) = c1*c1*c1;
    xi(21) =   a; eta(21) = 0.0; zeta(21) =  -a; weight(21) = c1*c1*c2;
    xi(22) =   a; eta(22) = 0.0; zeta(22) = 0.0; weight(22) = c1*c2*c2;
    xi(23) =   a; eta(23) = 0.0; zeta(23) =   a; weight(23) = c1*c1*c2;
    xi(24) =   a; eta(24) =   a; zeta(24) =  -a; weight(24) = c1*c1*c1;
    xi(25) =   a; eta(25) =   a; zeta(25) = 0.0; weight(25) = c1*c1*c2;
    xi(26) =   a; eta(26) =   a; zeta(26) =   a; weight(26) = c1*c1*c1;
}

void gauss_points_tetra1(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta){
    weight.Resize(1);
    xi.Resize(1); eta.Resize(1); zeta.Resize(1);
    weight(0) = 1.0/6.0;
    xi(0) = 1.0/4.0;
    eta(0) = 1.0/4.0;
    zeta(0) = 1.0/4.0;
}
void gauss_points_tetra4(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta){
    weight.Resize(4); xi.Resize(4); eta.Resize(4); zeta.Resize(4);
    double alpha = (5.0 - sqrt(5.0))/20.0;
    double beta = (5.0 + 3.0*sqrt(5.0))/20.0;
    weight(0) = 1.0/24.0; weight(1) = 1.0/24.0; weight(2) = 1.0/24.0; weight(3) = 1.0/24.0;
    xi(0) = alpha; eta(0) = alpha; zeta(0) = alpha;
    xi(1) = alpha; eta(1) = alpha; zeta(1) = beta;
    xi(2) = alpha; eta(2) = beta;  zeta(2) = alpha;
    xi(3) = beta;  eta(3) = alpha; zeta(3) = alpha;
}
void gauss_points_tetra5(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta){
    weight.Resize(5);xi.Resize(5); eta.Resize(5); zeta.Resize(5);
    xi(0) = 0.25; xi(1) = 0.50; xi(3) = 0.1666666666666667; xi(4) = 0.1666666666666667; xi(5) = 0.1666666666666667;
    eta(0) = 0.25; eta(1) = 0.1666666666666667; eta(2) = 0.1666666666666667; eta(3) = 0.1666666666666667; eta(4) = 0.50;
    zeta(0) = 0.25; zeta(1) = 0.1666666666666667; zeta(2) = 0.1666666666666667; zeta(3) = 0.50; zeta(4) = 0.1666666666666667;
    weight(0) = -0.8000000000000000/6.0; weight(1) = 0.45/6.0; weight(2) = 0.45/6.0; weight(3) = 0.45/6.0; weight(4) =0.45/6.0;
    
}
void gauss_points_tetra11(Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta){
    weight.Resize(11); xi.Resize(11); eta.Resize(11); zeta.Resize(11);
    xi(0) = 0.25; xi(1) = 0.7857142857142857; xi(2) = 0.0714285714285714; xi(3) = 0.0714285714285714; xi(4) = 0.0714285714285714; xi(5) = 0.1005964238332008; xi(6) = 0.3994035761667992; xi(7) = 0.3994035761667992; xi(8) = 0.3994035761667992; xi(9) = 0.1005964238332008; xi(10) = 0.1005964238332008;
    eta(0) = 0.25; eta(1) = 0.0714285714285714; eta(2) = 0.0714285714285714; eta(3) = 0.0714285714285714; eta(4) = 0.7857142857142857; eta(5) = 0.3994035761667992; eta(6) = 0.1005964238332008; eta(7) = 0.3994035761667992; eta(8) = 0.1005964238332008; eta(9) = 0.3994035761667992; eta(10) = 0.1005964238332008;
    zeta(0) = 0.25; zeta(1) = 0.0714285714285714; zeta(2) = 0.0714285714285714; zeta(3) = 0.7857142857142857; zeta(4) = 0.0714285714285714; zeta(5) = 0.3994035761667992; zeta(6) = 0.3994035761667992; zeta(7) = 0.1005964238332008; zeta(8) = 0.1005964238332008; zeta(9) = 0.1005964238332008; zeta(10) = 0.3994035761667992;
    weight(0) = -0.0789333333333333/6.0; weight(1) = 0.0457333333333333/6.0; weight(2) = 0.0457333333333333/6.0; weight(3) = 0.0457333333333333/6.0; weight(4) = 0.0457333333333333/6.0; weight(5) = 0.1493333333333333/6.0; weight(6) = 0.1493333333333333/6.0; weight(7) = 0.1493333333333333/6.0; weight(8) = 0.1493333333333333/6.0; weight(9) = 0.1493333333333333/6.0; weight(10) = 0.1493333333333333/6.0;
}

int pnpoly(int & nvert, Epetra_SerialDenseVector & vertx, Epetra_SerialDenseVector & verty, double & testx, double & testy){
    int i, j, c = 0;
    double a, b;
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
        a = (verty[j]-verty[i])/(vertx[j]-vertx[i]);
        b = verty[i]-a*vertx[i];
        if (testy - a*testx - b == 0){
            return 1;
        }
        if ( ((verty[i]>testy) != (verty[j]>testy)) &&
            (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ){
            c = !c;
        }
        
    }
    return c;
}
