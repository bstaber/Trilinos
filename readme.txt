Things that need to be improved:
1. The mesh class should be broken in multiple files.
2. The namespaces in tri3, quad4, tetra3, tetra10, hexa8, etc. in fepp.cpp are annoying and useless imo. Since we are only using Lagrange polynomials, the code could be generalized without having recourse to these namespaces.
3. In meshpp.cpp, the coordinates are sometimes divided by 1000.0 because I used a mesh in mm. But when I am given a mesh in meters, I need to modify a few lines: inefficient.
4. Gauss quadratures points in fepp.
5. Neumann type boundary conditions are not general enough. Should add a virtual class. 
