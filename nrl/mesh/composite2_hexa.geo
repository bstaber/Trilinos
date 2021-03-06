// Gmsh project created on Thu Mar 16 15:55:53 2017

//DefineConstant[ N = {4, Name "Number of layers"} ];

Mesh.RecombineAll=1;
Mesh.RecombinationAlgorithm=1;
Mesh.Algorithm = 8; 

N = 32;
lc = 0.5;
lcc = 0.5;

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

Point(17) = {8,14,0,lc};
Point(18) = {8,11,0,lc};
Point(19) = {8,25,0,lc};
Point(20) = {8,0,0,lc};
Point(21) = {0,18,0,lc};
Point(22) = {0,7,0,lc};
Point(23) = {8,18,0,lc};
Point(24) = {8,7,0,lc};
Point(25) = {16,18,0,lc};
Point(26) = {16,7,0,lc};
Point(27) = {16,25,0,lc};
Point(28) = {16,0,0,lc};

Point(29) = {42,14,0,lc};
Point(30) = {42,11,0,lc};
Point(31) = {42,25,0,lc};
Point(32) = {42,0,0,lc};
Point(33) = {42,18,0,lc};
Point(34) = {42,7,0,lc};
Point(35) = {34,18,0,lc};
Point(36) = {34,7,0,lc};
Point(37) = {34,25,0,lc};
Point(38) = {34,0,0,lc};
Point(39) = {50,18,0,lc};
Point(40) = {50,7,0,lc};//+
Line(1) = {7, 19};
//+
Line(2) = {19, 27};
//+
Line(3) = {27, 37};
//+
Line(4) = {37, 31};
//+
Line(5) = {31, 14};
//+
Line(6) = {14, 39};
//+
Line(7) = {39, 13};
//+
Line(8) = {13, 29};
//+
Line(9) = {29, 12};
//+
Line(10) = {10, 30};
//+
Line(11) = {30, 9};
//+
Line(12) = {9, 40};
//+
Line(13) = {40, 8};
//+
Line(14) = {8, 32};
//+
Line(15) = {32, 38};
//+
Line(16) = {38, 28};
//+
Line(17) = {28, 20};
//+
Line(18) = {20, 1};
//+
Line(19) = {1, 22};
//+
Line(20) = {22, 2};
//+
Line(21) = {2, 18};
//+
Line(22) = {18, 3};
//+
Line(23) = {5, 17};
//+
Line(24) = {17, 6};
//+
Line(25) = {6, 21};
//+
Line(26) = {21, 7};
//+
Line(27) = {17, 23};
//+
Line(28) = {23, 21};
//+
Line(29) = {23, 19};
//+
Line(30) = {23, 25};
//+
Line(31) = {25, 27};
//+
Line(32) = {37, 35};
//+
Line(33) = {35, 33};
//+
Line(34) = {33, 31};
//+
Line(35) = {33, 39};
//+
Line(36) = {33, 29};
//+
Line(37) = {18, 24};
//+
Line(38) = {24, 22};
//+
Line(39) = {24, 20};
//+
Line(40) = {24, 26};
//+
Line(41) = {26, 28};
//+
Line(42) = {25, 35};
//+
Line(43) = {26, 36};
//+
Line(44) = {36, 38};
//+
Line(45) = {36, 34};
//+
Line(46) = {34, 32};
//+
Line(47) = {34, 30};
//+
Line(48) = {34, 40};
//+
Circle(49) = {5, 4, 15};
//+
Circle(50) = {15, 4, 3};
//+
Circle(51) = {12, 11, 16};
//+
Circle(52) = {16, 11, 10};
//+
Line Loop(53) = {26, 1, -29, 28};
//+
Plane Surface(54) = {53};
//+
Line Loop(55) = {28, -25, -24, 27};
//+
Plane Surface(56) = {55};
//+
Line Loop(57) = {29, 2, -31, -30};
//+
Plane Surface(58) = {57};
//+
Line Loop(59) = {31, 3, 32, -42};
//+
Plane Surface(60) = {59};
//+
Line Loop(61) = {32, 33, 34, -4};
//+
Plane Surface(62) = {61};
//+
Line Loop(63) = {34, 5, 6, -35};
//+
Plane Surface(64) = {63};
//+
Line Loop(65) = {36, -8, -7, -35};
//+
Plane Surface(66) = {65};
//+
Line Loop(67) = {11, 12, -48, 47};
//+
Plane Surface(68) = {67};
//+
Line Loop(69) = {46, -14, -13, -48};
//+
Plane Surface(70) = {69};
//+
Line Loop(71) = {45, 46, 15, -44};
//+
Plane Surface(72) = {71};
//+
Line Loop(73) = {43, 44, 16, -41};
//+
Plane Surface(74) = {73};
//+
Line Loop(75) = {40, 41, 17, -39};
//+
Plane Surface(76) = {75};
//+
Line Loop(77) = {38, -19, -18, -39};
//+
Plane Surface(78) = {77};
//+
Line Loop(79) = {38, 20, 21, 37};
//+
Plane Surface(80) = {79};
//+
Line(81) = {25, 26};
//+
Line(82) = {36, 35};
//+
Line Loop(83) = {42, -82, -43, -81};
//+
Plane Surface(84) = {83};
//+
Line Loop(85) = {33, 36, 9, 51, 52, 10, -47, -45, 82};
//+
Plane Surface(86) = {85};
//+
Line Loop(87) = {27, 30, 81, -40, -37, 22, -50, -49, 23};
//+
Plane Surface(88) = {87};
//+
Extrude {0, 0, 4.10} {
  Surface{54, 58, 60, 62, 64, 66, 56, 88, 80, 78, 76, 84, 74, 72, 86, 68, 70};
  Layers{N};
  Recombine;
}
//+
Physical Surface(513) = {328, 350, 394, 416, 503, 189, 175, 145, 123, 101};
//+
Physical Surface(514) = {54, 110, 58, 132, 60, 154, 62, 176, 64, 198, 56, 242, 289, 88, 80, 311, 78, 333, 76, 355, 377, 84, 399, 74, 468, 86, 72, 421, 512, 70, 490, 68, 220, 66};
//+
Physical Volume(515) = {1, 2, 3, 4, 5, 6, 15, 12, 8, 7, 9, 10, 11, 13, 14, 17, 16};
