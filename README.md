# Three-dimensional Finite Element Method  

* Build on top of the [Trilinos project](https://github.com/trilinos/Trilinos)

* Uses: [GMSH](http://gmsh.info/), [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview), [Boost](https://www.boost.org/), [Epetra](https://trilinos.org/packages/epetra/), [Amesos](https://trilinos.org/packages/amesos/), [AztecOO](https://trilinos.org/packages/aztecoo/), [Stratimikos](https://trilinos.org/packages/stratimikos/), [Teuchos](https://trilinos.org/packages/teuchos/), [ML](https://trilinos.org/packages/ml/), amongst others.
* Handles linear and quadratic interpolations (tetra, hexas, triangle, quads) and various Gauss quadratures
* Linearized elasticity
* Compressible hyperelasticity
* Nearly incompressible hyperelasticity (three-field formulation with static condensation)
* Also solves Laplace's equation for the generation of fibers directions
* Essential and natural boundary conditions, deformation dependent natural conditions
* Stochastic boundary value problems with random field models
* Simulation of arterial walls with uncertainties (generation of Gaussian fields not included)
* [Documentation](https://bstaber.github.io/Trilinos/inherits.html)


