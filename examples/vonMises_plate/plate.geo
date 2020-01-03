// Gmsh project created on Mon May 22 16:45:07 2017

lc = 4.0;
Point(0) = {0,0,0,lc};
Point(1) = {4,0,0,lc};
Point(2) = {4,10,0,lc};
Point(3) = {0,10,0,lc};

Point(4) = {0,0,0.2,lc};
Point(5) = {4,0,0.2,lc};
Point(6) = {4,10,0.2,lc};
Point(7) = {0,10,0.2,lc};//+
Line(1) = {4, 5};
//+
Line(2) = {5, 6};
//+
Line(3) = {6, 7};
//+
Line(4) = {7, 4};
//+
Line(5) = {4, 0};
//+
Line(6) = {0, 1};
//+
Line(7) = {1, 5};
//+
Line(8) = {1, 2};
//+
Line(9) = {6, 2};
//+
Line(10) = {2, 3};
//+
Line(11) = {3, 7};
//+
Line(12) = {3, 0};
//+
Line Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {6, 8, 10, 12};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {8, -9, -2, -7};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {9, 10, 11, -3};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {12, -5, -4, -11};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {1, -7, -6, -5};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {1, 3, 2, 6, 5, 4};
//+
Volume(1) = {1};
//+
Physical Surface(1) = {4};
//+
Physical Surface(2) = {6};
//+
Physical Volume(3) = {1};
