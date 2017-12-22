Things that need to be improved:
1. The mesh class should be broken in multiple files: shape class, etc.
2. The namespaces in tri3, quad4, tetra3, tetra10, hexa8, etc. in fepp.cpp are annoying and useless imo. Since we are only using Lagrange polynomials, the code could be generalized without having recourse to these namespaces.
3. Units in meshpp.cpp.
4. Gauss quadratures points in fepp.
