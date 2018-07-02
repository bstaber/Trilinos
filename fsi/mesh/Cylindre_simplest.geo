cl__1 = 1;
//+
factor = 0.0254;
//+
Rmean = 1.485*factor;
//=+
ecyl = 0.0099*factor; //0.0099
//+
Rint = Rmean-ecyl/2; 
//+
Rext = Rmean+ecyl/2; 
//+
hdisque = 0.001; //0.001
//+
hinterieur = 9*factor;
hliquid = 0.5*hinterieur;
//+
ellipse_major = 0.25*factor;
ellipse_center = 0.25*factor;
//+
//+
hcapi_inf = 0.05*factor;
hcapi_mid = 0.1775*factor; //0.1567*factor;
hcapi_sup = 0.18745*factor;
//+
hvide = hinterieur - hliquid - hcapi_sup + hcapi_inf;
//+
distmaill = 1*factor; //0.035; //
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, Rext, 0, 1.0};
//+
Point(3) = {0, Rint, 0, 1.0};
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Duplicata { Point{2, 3}; }
}

//+
//+
//+
Line(1) = {1, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {1, 5};
//+
Line(4) = {5, 4};
//+
Circle(5) = {3, 1, 5};
//+
Circle(6) = {2, 1, 4};
//+
Line Loop(1) = {3, -5, -1};
//+
Surface(1) = {1};
//+
Line Loop(2) = {5, 4, -6, -2};
//+
Surface(2) = {2};
//+
Extrude {0, 0, hdisque} {
  Surface{1}; Surface{2}; 
}
//+
Extrude {0, 0, hliquid-hcapi_inf} {
  Surface{23}; Surface{45}; 
}
//+
Extrude {0, 0, hcapi_sup} {
  Surface{84}; 
}
//+
Translate {0, 0, hcapi_inf} {
  Duplicata { Point{25}; }
}
//+
Translate {0, 0, ellipse_center} {
  Duplicata { Point{25}; }
}
//+
Translate {0, distmaill, 0} {
  Duplicata { Point{57}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Duplicata { Point{58}; }
}
//+
Ellipse(101) = {56, 57, 59, 46};
//+
Ellipse(102) = {56, 57, 58, 44};
//+
Line(103) = {25, 56};
//+
Line Loop(3) = {47, 92, -101, -103};
//+
Surface(107) = {3};
//+
Line Loop(4) = {103, 102, -91, 49};
//+
Surface(108) = {4};
//+
Line Loop(5) = {101, -86, -102};
//+
Surface(109) = {5};
//+
Extrude {0, 0, hvide} {
  Surface{106}; 
}
//+q
Line(126) = {61, 60};
//+
Line(127) = {61, 62};
//+
Line Loop(6) = {127, -111, -126};
//+
Surface(132) = {6};
//+
Extrude {0, 0, hdisque} {
  Surface{132}; Surface{131}; 
}
//+
Surface Loop(1) = {108, 107, 109, 62, 93};
//+
Volume(9) = {1};
//+
Surface Loop(2) = {149, 140, 132, 148, 144};
//+
Volume(10) = {2};

//+
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Duplicata { Volume{7}; Volume{8}; Volume{10}; Volume{6}; Volume{5}; Volume{9}; Volume{3}; Volume{4}; Volume{1}; Volume{2}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Duplicata { Volume{172}; Volume{196}; Volume{227}; Volume{251}; Volume{282}; Volume{313}; Volume{337}; Volume{361}; Volume{392}; Volume{416}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Duplicata { Volume{438}; Volume{462}; Volume{493}; Volume{517}; Volume{548}; Volume{579}; Volume{603}; Volume{627}; Volume{658}; Volume{682}; }
}

//+
Physical Volume("solid", 1) = {196, 172, 227, 462, 438, 493, 7, 10, 8, 704, 759, 728, 251, 517, 6, 783, 416, 392, 682, 658, 1, 2, 924, 948, 548, 282, 5, 814, 627, 893, 4, 361};
//+
Physical Volume("fluid", 2) = {579, 845, 9, 313, 603, 337, 3, 869};
//+
Physical Surface("freesurf", 7) = {590, 856, 324, 109};
//+
Physical Surface("surfcoup", 9) = {351, 617, 883, 57, 293, 559, 825, 93, 338, 604, 870, 23};
//+
Physical Line("tripleline", 8) = {253, 519, 785, 86};


Nhdisk = 2;
Nhvide = 39; //13;
Nhcapi = 2; //;;
Nhfluid = 42; //7
Nepcyl = 2;

Nperimetre = 35; //8
Nperimetre_int = 15; //8
Nrayon = 15; //7
Nrayon_int = 11;
//+
Transfinite Line {535, 530, 269, 264, 116, 125, 117, 121} = Nhvide Using Progression 1;
//+
Transfinite Line {56, 78, 645, 614, 379, 348, 52, 74} = Nhfluid Using Progression 1;
//+
Transfinite Line {51} = 8 Using Progression 1;
//+
Transfinite Line {700, 669, 434, 12, 17, 39, 403, 13, 35, 480, 449, 138, 214, 183, 165, 143, 139, 161} = Nhdisk Using Progression 1;
//+
Transfinite Line {91, 100, 295, 300, 92, 96, 561, 566} = Nhcapi Using Progression 1;
//+
Transfinite Line {103} = 2 Using Progression 1;
//+
Transfinite Line {711, 737, 732, 707, 785, 787, 818, 816, 897, 872, 952, 927, 153, 135, 111, 113, 86, 88, 66, 48, 27, 9, 5, 6, 179, 205, 200, 175, 253, 255, 286, 284, 365, 340, 420, 395, 445, 471, 441, 466, 521, 519, 550, 552, 631, 606, 686, 661} = Nperimetre Using Progression 1;
//+
Transfinite Line {444, 440, 605, 660, 8, 3, 339, 394, 10, 1, 134, 127, 136, 126, 178, 174} = Nrayon Using Progression 1.25;
Transfinite Line {444, 440, 605, 660, 8, 3, 339, 394, 1, 134, 127, 126, 178, 174} = Nrayon Using Progression 0.8;
//+
Transfinite Line {588, 586, 102, 49, 101, 47, 322, 320} = Nrayon_int Using Progression 1;
