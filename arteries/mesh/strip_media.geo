// Gmsh project created on Mon May 22 16:45:07 2017

lc = 0.08;
Point(0) = {0,0,0,lc};
Point(1) = {3,0,0,lc};
Point(2) = {3,10,0,lc};
Point(3) = {0,10,0,lc};

Point(4) = {0,0,0.5,lc};
Point(5) = {3,0,0.5,lc};
Point(6) = {3,10,0.5,lc};
Point(7) = {0,10,0.5,lc};//+
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
Line Loop(13) = {4, 1, 2, 3};
//+
Plane Surface(14) = {13};
//+
Line Loop(15) = {3, -11, -10, -9};
//+
Plane Surface(16) = {15};
//+
Line Loop(17) = {8, -9, -2, -7};
//+
Plane Surface(18) = {17};
//+
Line Loop(19) = {6, 7, -1, 5};
//+
Plane Surface(20) = {19};
//+
Line Loop(21) = {4, 5, -12, 11};
//+
Plane Surface(22) = {21};
//+
Line Loop(23) = {8, 10, 12, 6};
//+
Plane Surface(24) = {23};
//+
Surface Loop(25) = {14, 22, 20, 24, 18, 16};
//+
Volume(26) = {25};
//+
Physical Surface(27) = {20};
//+
Physical Surface(1) = {16};
//+
Physical Surface(29) = {22};
//+
Physical Surface(30) = {18};
//+
Physical Surface(31) = {14};
//+
Physical Surface(32) = {24};
//+
Physical Volume(33) = {26};
//+
//Transfinite Line {11, 9, 5, 7} = 20 Using Progression 1;
//+
//Transfinite Line {3, 10, 1, 6} = 60 Using Progression 1;
//+
//Transfinite Line {12, 4, 8, 2} = 100 Using Progression 1;
