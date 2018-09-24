// Gmsh project created on Sun Sep 23 18:09:47 2018
SetFactory("OpenCASCADE");

lc1 = 0.25;
lc2 = 0.05;

H = 0.1;

a = 2;
b = 5;
c = 10;
d = 3;

Point(1) = {0, 0, 0, lc1};
Point(2) = {b, 0, 0, lc1};
Point(3) = {0, a, 0, lc1};
Point(4) = {b, a, 0, lc1};

Point(5) = {0, 0, c, lc1};
Point(6) = {b, 0, c, lc1};
Point(7) = {0, a, c, lc1};
Point(8) = {b, a, c, lc1};

Point(9)  = {d, 0, b, lc2};
Point(10) = {b, 0, b-H, lc2};
Point(11) = {d, a, b, lc2};
Point(12) = {b, a, b-H, lc2};

Point(13) = {b, 0, b+H, lc2};
Point(14) = {b, a, b+H, lc2};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Line(5) = {4, 12};
//+
Line(6) = {12, 10};
//+
Line(7) = {10, 2};
//+
Line(8) = {10, 9};
//+
Line(9) = {9, 11};
//+
Line(10) = {11, 12};
//+
Line(11) = {9, 13};
//+
Line(12) = {11, 14};
//+
Line(13) = {14, 13};
//+
Line(14) = {14, 8};
//+
Line(15) = {8, 6};
//+
Line(16) = {6, 13};
//+
Line(17) = {8, 7};
//+
Line(18) = {7, 5};
//+
Line(19) = {5, 6};
//+
Line(20) = {5, 1};
//+
Line(21) = {3, 7};
//+
Line Loop(1) = {21, -17, -14, -12, 10, -5, 3};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {19, 16, -11, -8, 7, -1, -20};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {19, -15, 17, 18};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {21, 18, 20, -4};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {1, 2, 3, 4};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {2, 5, 6, 7};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {16, -13, 14, 15};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {12, 13, -11, 9};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {8, 9, 10, 6};
//+
Plane Surface(9) = {9};
//+
Surface Loop(1) = {2, 3, 7, 8, 1, 4, 5, 6, 9};
//+
Volume(1) = {1};
//+
Physical Surface(1) = {3};
//+
Physical Volume(2) = {1};
