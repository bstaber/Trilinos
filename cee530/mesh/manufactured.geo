// Gmsh project created on Thu Mar 16 15:55:53 2017

DefineConstant[ N = {4, Name "Number of layers"} ];

Mesh.RecombineAll=1;
Mesh.RecombinationAlgorithm=1;
Mesh.Algorithm = 8; 

Mesh.CharacteristicLengthMax = 4.1;

//lc = 8.0;
//lcc = 8.0;

lc = 1000;
lcc = 1000;

Point(1) = {0,0,0,lc};
Point(2) = {0,11,0,lc};
Point(3) = {10.5,11,0,lcc};
Point(4) = {10.5,12.5,0,lcc};
Point(5) = {10.5,14,0,lcc};
Point(6) = {0,14,0,lc};
Point(7) = {0,25,0,lc};
Point(15) = {12,12.5,0,lcc};

Point(8) = {50,0,0,lc};
Point(9) = {50,11,0,lc};
Point(10) = {39.5,11,0,lcc};
Point(11) = {39.5,12.5,0,lcc};
Point(12) = {39.5,14,0,lcc};
Point(13) = {50,14,0,lc};
Point(14) = {50,25,0,lc};
Point(16) = {38,12.5,0,lcc};
//+
Line(1) = {1, 8};
//+
Line(2) = {8, 9};
//+
Line(3) = {9, 10};
//+
Line(4) = {1, 2};
//+
Line(5) = {2, 3};
//+
Line(6) = {7, 6};
//+
Line(7) = {6, 5};
//+
Line(8) = {7, 14};
//+
Line(9) = {14, 13};
//+
Line(10) = {13, 12};
//+
Circle(11) = {5, 4, 15};
//+
Circle(12) = {15, 4, 3};
//+
Circle(13) = {12, 11, 16};
//+
Circle(14) = {16, 11, 10};
//+
Line Loop(15) = {8, 9, 10, 13, 14, -3, -2, -1, 4, 5, -12, -11, -7, -6};
//+
Plane Surface(16) = {-15};
//+
out[] = Extrude {0, 0, 4.10} {
    Surface{16};
    Layers{N};
    Using Index[0];
    Recombine;
  };

Physical Surface(92) = {35, 55, 39, 51, 43, 47, 88, 16, 75, 71, 79, 67, 83, 63};
Physical Volume(101) = {out[1]};
