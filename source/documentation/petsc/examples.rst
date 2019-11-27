.. _petscExamples:

Examples
--------

Linear problems
~~~~~~~~~~~~~~~

+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| Filename                                                                                                                                                       | Comments (preconditioners, numerical schemes)                              |
+================================================================================================================================================================+============================================================================+
| `diffusion-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-2d-PETSc.edp>`__                                     | Distributed LU/Cholesky, domain decomposition and multigrid methods        |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `diffusion-2d-PETSc-complex.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-2d-PETSc-complex.edp>`__                     |                                                                            |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `heat-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/heat-2d-PETSc.edp>`__                                               | Transient diffusion equation, same as above                                |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `diffusion-periodic-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-periodic-2d-PETSc.edp>`__                   | Periodic boundary conditions, multigrid methods                            |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `diffusion-periodic-balanced-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-periodic-balanced-2d-PETSc.edp>`__ | Better load balancing than above example                                   |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `diffusion-substructuring-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-substructuring-2d-PETSc.edp>`__       | Balancing Domain Decomposition with Constraints                            |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `diffusion-3d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-3d-PETSc.edp>`__                                     | Three-dimensional problem, domain decomposition and multigrid methods      |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `diffusion-mg-3d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-mg-3d-PETSc.edp>`__                               | Geometric multigrid methods                                                |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `helmholtz-2d-PETSc-complex.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-2d-PETSc-complex.edp>`__                     | Domain decomposition methods with optimized boundary conditions            |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `laplace-RT-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-RT-2d-PETSc.edp>`__                                   | Vectorial two-dimensional problem with a block preconditioner (fieldsplit) |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `laplace-adapt-3d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-adapt-3d-PETSc.edp>`__                             | Three-dimensional problem with *h* adaptivity, multigrid methods           |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `laplace-lagrange-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-lagrange-PETSc.edp>`__                             | Laplace equation with constraints and a block preconditioner (fieldsplit)  |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `elasticity-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/elasticity-2d-PETSc.edp>`__                                   | Vectorial problem, domain decomposition and multigrid methods              |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `elasticity-3d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/elasticity-3d-PETSc.edp>`__                                   |                                                                            |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `stokes-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-2d-PETSc.edp>`__                                           | Distributed LU/Cholesky                                                    |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `stokes-3d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-3d-PETSc.edp>`__                                           |                                                                            |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `stokes-block-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-block-2d-PETSc.edp>`__                               | Stokes equation defined as a block system with four matrices (fieldsplit)  |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `stokes-fieldsplit-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-fieldsplit-2d-PETSc.edp>`__                     | Block preconditioner (fieldsplit)                                          |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `stokes-fieldsplit-3d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-fieldsplit-3d-PETSc.edp>`__                     |                                                                            |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `maxwell-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/maxwell-2d-PETSc.edp>`__                                         | Direct LU/Cholesky                                                         |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+
| `maxwell-3d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/maxwell-3d-PETSc.edp>`__                                         | Multigrid method                                                           |
+----------------------------------------------------------------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------+

Nonlinear problems
~~~~~~~~~~~~~~~~~~

+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| Filename                                                                                                                                                           | Comments (preconditioners, numerical schemes)               |
+====================================================================================================================================================================+=============================================================+
| `bratu-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/bratu-2d-PETSc.edp>`__                                                 |                                                             |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `Newton-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/Newton-2d-PETSc.edp>`__                                               |                                                             |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `Newton-adaptmesh-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/Newton-adaptmesh-2d-PETSc.edp>`__                           | Newton method and *h* adaptivity                            |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `Newton-vi-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/Newton-vi-2d-PETSc.edp>`__                                         | Newton method and a variational inequality                  |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `Newton-vi-adaptmesh-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/Newton-vi-adaptmesh-2d-PETSc.edp>`__                     | Newton method, *h* adaptivity, and a variational inequality |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `elasticity-SNES-3d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/elasticity-SNES-3d-PETSc.edp>`__                             | Linear elasiticty with a Newton method                      |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `neo-Hookean-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/neo-Hookean-2d-PETSc.edp>`__                                     | Nonlinear elasticity                                        |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `natural-convection-fieldsplit-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/natural-convection-fieldsplit-2d-PETSc.edp>`__ | Newton method and *h* adaptivity                            |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `vi-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/vi-2d-PETSc.edp>`__                                                       | Variational inequalities                                    |
+--------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+

Time steppers and optimizers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| Filename                                                                                                                                       | Comments (preconditioners, numerical schemes) |
+================================================================================================================================================+===============================================+
| `advection-TS-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/advection-TS-2d-PETSc.edp>`__               | Implicit and explicit schemes                 |
+------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `heat-TS-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/heat-TS-2d-PETSc.edp>`__                         |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `heat-TS-RHS-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/heat-TS-RHS-2d-PETSc.edp>`__                 |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `minimal-surface-Tao-2d-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/minimal-surface-Tao-2d-PETSc.edp>`__ | Minimal surface problem                       |
+------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `orego-Tao-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/orego-Tao-PETSc.edp>`__                           |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `toy-Tao-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/toy-Tao-PETSc.edp>`__                               |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+

Eigenvalue problems
~~~~~~~~~~~~~~~~~~~

+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| Filename                                                                                                                                                               | Comments (preconditioners, numerical schemes) |
+========================================================================================================================================================================+===============================================+
| `laplace-2d-SLEPc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-2d-SLEPc.edp>`__                                                 |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `laplace-spherical-harmonics-2d-SLEPc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-spherical-harmonics-2d-SLEPc.edp>`__         |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `laplace-torus-2d-SLEPc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-torus-2d-SLEPc.edp>`__                                     |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `schrodinger-axial-well-2d-SLEPc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-axial-well-2d-SLEPc.edp>`__                   |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `schrodinger-harmonic-oscillator-1d-SLEPc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-harmonic-oscillator-1d-SLEPc.edp>`__ |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `schrodinger-harmonic-oscillator-2d-SLEPc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-harmonic-oscillator-2d-SLEPc.edp>`__ |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `schrodinger-square-well-1d-SLEPc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-square-well-1d-SLEPc.edp>`__                 |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+
| `laplace-2d-SLEPc-complex.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-2d-SLEPc-complex.edp>`__                                 |                                               |
+------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------+

Miscellaneous
~~~~~~~~~~~~~

+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| Filename                                                                                                                           | Comments (preconditioners, numerical schemes)               |
+====================================================================================================================================+=============================================================+
| `transpose-solve-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/transpose-solve-PETSc.edp>`__   | Solving a transposed system                                 |
+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `Schur-complement-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/Schur-complement-PETSc.edp>`__ | Computing an exact Schur complement                         |
+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `block-PETSc.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/block-PETSc.edp>`__                       |                                                             |
+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `heat-io-2d.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/heat-io-2d.edp>`__                         | Automatic ParaView animation output                         |
+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `stokes-io-3d.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-io-3d.edp>`__                     |                                                             |
+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `buildRecursive.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/buildRecursive.edp>`__                 | Recursive mesh partitioning (for geometric multigrid)       |
+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+
| `withPartitioning.edp <https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/withPartitioning.edp>`__             | Connectivity construction with a user-supplied partitioning |
+------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------+

Reproducible science
~~~~~~~~~~~~~~~~~~~~

+---------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| Article                                                                                                                                                       | Source code                                                     |
+===============================================================================================================================================================+=================================================================+
| `Augmented Lagrangian Preconditioner for Large-Scale Hydrodynamic Stability Analysis <https://www.sciencedirect.com/science/article/pii/S0045782519301914>`__ | `GitHub repository <https://github.com/prj-/moulin2019al>`__    |
+---------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| `A Multilevel Schwarz Preconditioner Based on a Hierarchy of Robust Coarse Spaces <https://hal.archives-ouvertes.fr/hal-02151184/document>`__                 | `GitHub repository <https://github.com/prj-/aldaas2019multi>`__ |
+---------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+

