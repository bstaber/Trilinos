/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef TENSOR_CALCULUS_HPP
#define TENSOR_CALCULUS_HPP

#include "meshpp.hpp"

void tensor_product(double scalarAB, Epetra_SerialDenseVector & A, Epetra_SerialDenseVector & B, Epetra_SerialDenseMatrix & AoB, double scalarThis){

    AoB.Multiply('N','T',scalarAB,A,B,scalarThis);

}

void sym_tensor_product(double scalarAB, Epetra_SerialDenseVector & A, Epetra_SerialDenseVector & B, Epetra_SerialDenseMatrix & AoB, double scalarThis){

    AoB(0,0) = AoB(0,0)*scalarThis + A(0)*B(0)*scalarAB; AoB(0,1) = AoB(0,1)*scalarThis + A(5)*B(5)*scalarAB; AoB(0,2) = AoB(0,2)*scalarThis + A(4)*B(4)*scalarAB; AoB(0,3) = AoB(0,3)*scalarThis + scalarAB*((A(4)*B(5))/2.0 + (A(5)*B(4))/2.0); AoB(0,4) = AoB(0,4)*scalarThis + scalarAB*((A(0)*B(4))/2.0 + (A(4)*B(0))/2.0); AoB(0,5) = AoB(0,5)*scalarThis + scalarAB*((A(0)*B(5))/2.0 + (A(5)*B(0))/2.0);
    AoB(1,0) = AoB(1,0)*scalarThis + A(5)*B(5)*scalarAB; AoB(1,1) = AoB(1,1)*scalarThis + A(1)*B(1)*scalarAB; AoB(1,2) = AoB(1,2)*scalarThis + A(3)*B(3)*scalarAB; AoB(1,3) = AoB(1,3)*scalarThis + scalarAB*((A(1)*B(3))/2.0 + (A(3)*B(1))/2.0); AoB(1,4) = AoB(1,4)*scalarThis + scalarAB*((A(3)*B(5))/2.0 + (A(5)*B(3))/2.0); AoB(1,5) = AoB(1,5)*scalarThis + scalarAB*((A(1)*B(5))/2.0 + (A(5)*B(1))/2.0);
    AoB(2,0) = AoB(2,0)*scalarThis + A(4)*B(4)*scalarAB; AoB(2,1) = AoB(2,1)*scalarThis + A(3)*B(3)*scalarAB; AoB(2,2) = AoB(2,2)*scalarThis + A(2)*B(2)*scalarAB; AoB(2,3) = AoB(2,3)*scalarThis + scalarAB*((A(2)*B(3))/2.0 + (A(3)*B(2))/2.0); AoB(2,4) = AoB(2,4)*scalarThis + scalarAB*((A(2)*B(4))/2.0 + (A(4)*B(2))/2.0); AoB(2,5) = AoB(2,5)*scalarThis + scalarAB*((A(3)*B(4))/2.0 + (A(4)*B(3))/2.0);
    AoB(3,0) = AoB(3,0)*scalarThis + A(5)*B(4)*scalarAB; AoB(3,1) = AoB(3,1)*scalarThis + A(1)*B(3)*scalarAB; AoB(3,2) = AoB(3,2)*scalarThis + A(3)*B(2)*scalarAB; AoB(3,3) = AoB(3,3)*scalarThis + scalarAB*((A(1)*B(2))/2.0 + (A(3)*B(3))/2.0); AoB(3,4) = AoB(3,4)*scalarThis + scalarAB*((A(3)*B(4))/2.0 + (A(5)*B(2))/2.0); AoB(3,5) = AoB(3,5)*scalarThis + scalarAB*((A(1)*B(4))/2.0 + (A(5)*B(3))/2.0);
    AoB(4,0) = AoB(4,0)*scalarThis + A(0)*B(4)*scalarAB; AoB(4,1) = AoB(4,1)*scalarThis + A(5)*B(3)*scalarAB; AoB(4,2) = AoB(4,2)*scalarThis + A(4)*B(2)*scalarAB; AoB(4,3) = AoB(4,3)*scalarThis + scalarAB*((A(4)*B(3))/2.0 + (A(5)*B(2))/2.0); AoB(4,4) = AoB(4,4)*scalarThis + scalarAB*((A(0)*B(2))/2.0 + (A(4)*B(4))/2.0); AoB(4,5) = AoB(4,5)*scalarThis + scalarAB*((A(0)*B(3))/2.0 + (A(5)*B(4))/2.0);
    AoB(5,0) = AoB(5,0)*scalarThis + A(0)*B(5)*scalarAB; AoB(5,1) = AoB(5,1)*scalarThis + A(5)*B(1)*scalarAB; AoB(5,2) = AoB(5,2)*scalarThis + A(4)*B(3)*scalarAB; AoB(5,3) = AoB(5,3)*scalarThis + scalarAB*((A(4)*B(1))/2.0 + (A(5)*B(3))/2.0); AoB(5,4) = AoB(5,4)*scalarThis + scalarAB*((A(0)*B(3))/2.0 + (A(4)*B(5))/2.0); AoB(5,5) = AoB(5,5)*scalarThis + scalarAB*((A(0)*B(1))/2.0 + (A(5)*B(5))/2.0);

}

#endif
