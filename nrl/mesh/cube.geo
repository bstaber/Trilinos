// Gmsh project created on Tue Mar 21 23:12:20 2017

lc = 0.05;

Point(0) = {0,0,0,lc};
Point(1) = {1,0,0,lc};
Point(2) = {1,1,0,lc};
Point(3) = {0,1,0,lc};//+
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
Transfinite Line {1, 4, 3, 2} = 10 Using Progression 1;
//+
//+
Recombine Surface {6};
//+
Transfinite Surface {6};
//+
Extrude {0, 0, 1} {
  Surface{6};
  Layers{10};
  Recombine;
}
