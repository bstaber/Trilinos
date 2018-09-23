// Gmsh project created on Sun Sep 23 18:09:47 2018
SetFactory("OpenCASCADE");

lc1 = 0.25;
lc2 = 0.05;

H = 0.001;

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

Point(9)  = {d, 0, b-H, lc2};
Point(10) = {b, 0, b-H, lc2};
Point(11) = {d, a, b-H, lc2};
Point(12) = {b, a, b-H, lc2};

Point(13) = {d, 0, b+H, lc2};
Point(14) = {b, 0, b+H, lc2};
Point(15) = {d, a, b+H, lc2};
Point(16) = {b, a, b+H, lc2};

//+
Line(1) = {2, 4};
//+
Line(2) = {2, 1};
//+
Line(3) = {1, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {4, 12};
//+
Line(6) = {12, 10};
//+
Line(7) = {10, 2};
//+
Line(8) = {10, 9};
//+
Line(9) = {9, 13};
//+
Line(10) = {13, 14};
//+
Line(11) = {14, 16};
//+
Line(12) = {16, 15};
//+
Line(13) = {15, 11};
//+
Line(14) = {11, 12};
//+
Line(15) = {9, 11};
//+
Line(16) = {13, 15};
//+
Line(17) = {16, 8};
//+
Line(18) = {8, 6};
//+
Line(19) = {6, 14};
//+
Line(20) = {8, 7};
//+
Line(21) = {7, 5};
//+
Line(22) = {5, 6};
//+
Line(23) = {7, 3};
//+
Line(24) = {5, 1};
//+
Line Loop(1) = {10, 11, 12, -16};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {16, 13, -15, 9};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {8, 15, 14, 6};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {17, 18, 19, 11};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {6, 7, 1, 5};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {1, -4, -3, -2};
//+
Plane Surface(6) = {6};
//+
Line Loop(7) = {23, -3, -24, -21};
//+
Plane Surface(7) = {7};
//+
Line Loop(8) = {22, -18, 20, 21};
//+
Plane Surface(8) = {8};
//+
Line Loop(9) = {2, -24, 22, 19, -10, -9, -8, 7};
//+
Plane Surface(9) = {9};
//+
Line Loop(10) = {23, 4, 5, -14, -13, -12, 17, 20};
//+
Plane Surface(10) = {10};
//+
Surface Loop(1) = {9, 6, 5, 3, 2, 1, 4, 10, 7, 8};
//+
Volume(1) = {1};
//+
Physical Surface("TOP") = {8};
//+
Physical Volume("VOL") = {1};
