// Gmsh project created on Tue Mar 21 23:12:20 2017
Mesh.RandomFactor=1.e-1;
lc = 5;
length = 10000;
height = 10000;

Point(0) = {0,0,0};
Point(1) = {length,0,0};
Point(2) = {length,0,length};
Point(3) = {0,0,length};//+
Line(1) = {0, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 0};
//+
Line Loop(5) = {4, 1, 2, 3};
//+
Plane Surface(6) = {5};
//+
Transfinite Line {1, 4, 3, 2} = lc Using Progression 1;
//+
//+
Recombine Surface {6};
//+
Transfinite Surface {6};
//+
Extrude {0, height, 0} {
  Surface{6};
  Layers{lc};
  Recombine;
}
//+
Physical Surface(1) = {19, 23, 28, 6, 15, 27};
//+
Physical Volume(2) = {1};
