# Three-dimensional Finite Element Method  

* Built on top of the [Trilinos project](https://github.com/trilinos/Trilinos)

* Uses: [GMSH](http://gmsh.info/), [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview), [Boost](https://www.boost.org/), [Epetra](https://trilinos.org/packages/epetra/), [Amesos](https://trilinos.org/packages/amesos/), [AztecOO](https://trilinos.org/packages/aztecoo/), [Stratimikos](https://trilinos.org/packages/stratimikos/), [Teuchos](https://trilinos.org/packages/teuchos/), [ML](https://trilinos.org/packages/ml/), amongst others.
* Handles **linear and quadratic interpolations** (tetra, hexas, triangle, quads) and various Gauss quadratures
* Most of the problems (arteries, asme, cee530, nrl) are **stochastic boundary value problems** where the elasticity tensor or strain energy function are **modeled as random fields**
* Deterministic problems are also considered in some subfolders
* The c++ classes defining the boundary value problems (i.e., **linearizedElasticity**, **compressibleHyperelasticity** and **nearlyIncompressibleHyperelasticity**) are pure virtual classes, **independent of the constitutive equations** (random or deterministic elasticity tensor, random or deterministic strain energy function), boundary conditions, etc.
	* __linearizedElasticity__: see folders asme, cee530
		* used in asmeSBVP, ceeSBVP, linearPatchTest and manufactured
	* __compressibleHyperelasticity__: see folder nrl
		* used in compressibleHyperelasticity_linearPatchTest, manufacturedSolution, rubberblock, tiMooney, tiMooneyRandomField
	* __nearlyIncompressibleHyperelasticity__ (three-field formulation with static condensation): see folder arteries
		* used in dirichletInletOutlet_PolyconvexHGO, dirichletStripElongation_StochasticPolyconvexHGO, neumannInnerSurface_PolyconvexHGO, neumannInnerSurface_StochasticPolyconvexHGO
* Also solves Laplace's equation for the generation of fibers directions in the folder arteries
* Handles essential and natural boundary conditions, deformation dependent natural conditions, with or without forcing terms.
* Have a look at the [class diagram](https://bstaber.github.io/Trilinos/inherits.html).


