Things that need to be improved:
1. The mesh class should be broken in multiple files: shape class, etc.
2. The namespaces in tri3, quad4, tetra3, tetra10, hexa8, etc. in fepp.cpp are annoying and useless imo. Since we are only using Lagrange polynomials, the code could be generalized without having recourse to these namespaces.
3. Units in meshpp.cpp.
4. Gauss quadratures points in fepp.
5. Need to add the line:
        <Parameter name="tol"               type="double" value="1e-8"/>
to all xml files with a Newton parameter list.

/!\ IMPORTANT /!\
6. Need to check the influence of AZ_tol (set to 1e-6 in the draft): try out 1e-8 and 1-10 for instance.
7. There was a tiny mistake in the tangent operator of gmrf_neumann and strips: not a big deal but it has been corrected.
