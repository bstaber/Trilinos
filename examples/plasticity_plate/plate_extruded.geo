// Gmsh project created on Mon May 22 16:45:07 2017

Mesh.RecombineAll=1;
Mesh.RecombinationAlgorithm=1;
Mesh.Algorithm = 8; 
Mesh.SecondOrderIncomplete = 1;

lc = 4.0;
Point(0) = {0,0,0,lc};
Point(1) = {4,0,0,lc};
Point(2) = {4,10,0,lc};
Point(3) = {0,10,0,lc};

Line(1) = {0, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 0};

//+
Transfinite Line {1} = 7 Using Progression 1;
//+
Transfinite Line {4} = 16 Using Progression 1;
//+
Transfinite Line {3} = 7 Using Progression 1;
//+
Transfinite Line {2} = 16 Using Progression 1;
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
out[] = Extrude {0, 0, 0.2} {
    Surface{1};
    Layers{1};
    Using Index[0];
    Recombine;
  };

Physical Surface(102) = {26, 1, 13, 17, 21, 25};
Physical Volume(101) = {out[1]};
