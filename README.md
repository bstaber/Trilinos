
## TrilinosUQComp: Uncertainty quantification in computational mechanics

### Packages and libraries

* Built on top of Sandia's Trilinos Project [[Website](http://trilinos.org/) |
[Documentation](http://trilinos.org/about/documentation/) |
[Mailing List](https://trilinos.org/mailman/listinfo/trilinos-users) |
[Packages](http://trilinos.org/packages/) |
[GitHub](https://github.com/trilinos/Trilinos)] and uses, for instance, the packages:
	* [Epetra](https://trilinos.org/packages/epetra/) for linear algebra,
	* [Amesos](https://trilinos.org/packages/amesos/), [AztecOO](https://trilinos.org/packages/aztecoo/), [Stratimikos](https://trilinos.org/packages/stratimikos/) for direct and iterative linear solvers,
	* [Teuchos](https://trilinos.org/packages/teuchos/) for parsing XML lists of parameters,
	* [ML](https://trilinos.org/packages/ml/), [IfPack](https://trilinos.org/packages/ifpack/) for preconditioners.

* Additional required libraries: [GMSH](http://gmsh.info/), [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview), [Boost](https://www.boost.org/).

### Documentation

* Have a look at the [class diagram](https://bstaber.github.io/Trilinos/inherits.html).

<!--
### Description

* Handles **linear and quadratic interpolations** (tetras, hexas, triangles, quads) and various Gauss quadratures.
* Handles essential and natural boundary conditions, deformation dependent natural conditions, with or without forcing terms.
* Most of the problems (arteries, asme, cee530, nrl) are **stochastic boundary value problems** where the elasticity tensor or strain energy function are **modeled as random fields**.
* Deterministic problems are also considered in some subfolders.
* The c++ classes defining the boundary value problems (i.e., **linearizedElasticity**, **compressibleHyperelasticity** and **nearlyIncompressibleHyperelasticity**) are pure virtual classes, **independent of the constitutive equations** (random or deterministic elasticity tensor, random or deterministic strain energy function), boundary conditions, etc.
	* __linearizedElasticity__: see folders asme, cee530
		* used in asmeSBVP, ceeSBVP, linearPatchTest and manufactured
	* __compressibleHyperelasticity__: see folder nrl
		* used in compressibleHyperelasticity_linearPatchTest, manufacturedSolution, rubberblock, tiMooney, tiMooneyRandomField
	* __nearlyIncompressibleHyperelasticity__ (three-field formulation with static condensation): see folder arteries
		* used in dirichletInletOutlet_PolyconvexHGO, dirichletStripElongation_StochasticPolyconvexHGO, neumannInnerSurface_PolyconvexHGO, neumannInnerSurface_StochasticPolyconvexHGO
* Also solves Laplace's equation for the generation of fibers directions in the folder arteries.

### Example: solving a boundary value problem in linearized elasticity

* The pure virtual class **linearizedElasticity** has six virtual methods, namely:
	* **get_neumannBc**: returns your natural boundary condition while assembling the force vector
	* **get_forcing**: returns your forcing vector while assembling
	* **get_elasticity_tensor**: returns your elasticity tensor
	* **get_elasticity_tensor_for_recovery**: returns the elasticity tensor for post-processing the solution, useful for stochastic boundary value problems
	* **setup_dirichlet_conditions**: builds the essential boundary conditions that you want to apply
	* **apply_dirichlet_conditions**: apply your essential boundary conditions

* You need to construct your own class object inheriting from **linearizedElasticity** and which contains at least the six virtual methods.

* Examples can be found in e.g. **src/asmeSBVP.hpp**, **cee530/sbvp/ceeSBVP.hpp**, **cee530/linearPatchTest/linearPatchTest.hpp** and **cee530/manufactured/manufactured.hpp**
-->
### Warning

* Line 167 in meshpp.cpp needs to be removed/changed.
