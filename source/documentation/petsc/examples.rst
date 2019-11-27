.. _petscExamples:

Examples
--------

Linear problems
~~~~~~~~~~~~~~~

+-----------------------------------------------+----------------------+
| Filename                                      | Comments             |
|                                               | (preconditioners,    |
|                                               | numerical schemes)   |
+===============================================+======================+
| `diffusion-2d-PETSc.edp <https:/              | Distributed          |
| /github.com/FreeFem/FreeFem-sources/tree/deve | LU/Cholesky, domain  |
| lop/examples/hpddm/diffusion-2d-PETSc.edp>`__ | decomposition and    |
|                                               | multigrid methods    |
+-----------------------------------------------+----------------------+
| `di                                           |                      |
| ffusion-2d-PETSc-complex.edp <https://github. |                      |
| com/FreeFem/FreeFem-sources/tree/develop/exam |                      |
| ples/hpddm/diffusion-2d-PETSc-complex.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `heat-2d-PETSc.edp <ht                        | Transient diffusion  |
| tps://github.com/FreeFem/FreeFem-sources/tree | equation, same as    |
| /develop/examples/hpddm/heat-2d-PETSc.edp>`__ | above                |
+-----------------------------------------------+----------------------+
| `diff                                         | Periodic boundary    |
| usion-periodic-2d-PETSc.edp <https://github.c | conditions,          |
| om/FreeFem/FreeFem-sources/tree/develop/examp | multigrid methods    |
| les/hpddm/diffusion-periodic-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `diffusion-periodic-bal                       | Better load          |
| anced-2d-PETSc.edp <https://github.com/FreeFe | balancing than above |
| m/FreeFem-sources/tree/develop/examples/hpddm | example              |
| /diffusion-periodic-balanced-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `diffusion-substr                             | Balancing Domain     |
| ucturing-2d-PETSc.edp <https://github.com/Fre | Decomposition with   |
| eFem/FreeFem-sources/tree/develop/examples/hp | Constraints          |
| ddm/diffusion-substructuring-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `diffusion-3d-PETSc.edp <https:/              | Three-dimensional    |
| /github.com/FreeFem/FreeFem-sources/tree/deve | problem, domain      |
| lop/examples/hpddm/diffusion-3d-PETSc.edp>`__ | decomposition and    |
|                                               | multigrid methods    |
+-----------------------------------------------+----------------------+
| `diffusion-mg-3d-PETSc.edp <https://gi        | Geometric multigrid  |
| thub.com/FreeFem/FreeFem-sources/tree/develop | methods              |
| /examples/hpddm/diffusion-mg-3d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `he                                           | Domain decomposition |
| lmholtz-2d-PETSc-complex.edp <https://github. | methods with         |
| com/FreeFem/FreeFem-sources/tree/develop/exam | optimized boundary   |
| ples/hpddm/helmholtz-2d-PETSc-complex.edp>`__ | conditions           |
+-----------------------------------------------+----------------------+
| `laplace-RT-2d-PETSc.edp <https://            | Vectorial            |
| github.com/FreeFem/FreeFem-sources/tree/devel | two-dimensional      |
| op/examples/hpddm/laplace-RT-2d-PETSc.edp>`__ | problem with a block |
|                                               | preconditioner       |
|                                               | (fieldsplit)         |
+-----------------------------------------------+----------------------+
| `laplace-adapt-3d-PETSc.edp <https://git      | Three-dimensional    |
| hub.com/FreeFem/FreeFem-sources/tree/develop/ | problem with *h*     |
| examples/hpddm/laplace-adapt-3d-PETSc.edp>`__ | adaptivity,          |
|                                               | multigrid methods    |
+-----------------------------------------------+----------------------+
| `laplace-lagrange-PETSc.edp <https://git      | Laplace equation     |
| hub.com/FreeFem/FreeFem-sources/tree/develop/ | with constraints and |
| examples/hpddm/laplace-lagrange-PETSc.edp>`__ | a block              |
|                                               | preconditioner       |
|                                               | (fieldsplit)         |
+-----------------------------------------------+----------------------+
| `elasticity-2d-PETSc.edp <https://            | Vectorial problem,   |
| github.com/FreeFem/FreeFem-sources/tree/devel | domain decomposition |
| op/examples/hpddm/elasticity-2d-PETSc.edp>`__ | and multigrid        |
|                                               | methods              |
+-----------------------------------------------+----------------------+
| `elasticity-3d-PETSc.edp <https://            |                      |
| github.com/FreeFem/FreeFem-sources/tree/devel |                      |
| op/examples/hpddm/elasticity-3d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `stokes-2d-PETSc.edp <http                    | Distributed          |
| s://github.com/FreeFem/FreeFem-sources/tree/d | LU/Cholesky          |
| evelop/examples/hpddm/stokes-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `stokes-3d-PETSc.edp <http                    |                      |
| s://github.com/FreeFem/FreeFem-sources/tree/d |                      |
| evelop/examples/hpddm/stokes-3d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `stokes-block-2d-PETSc.edp <https://gi        | Stokes equation      |
| thub.com/FreeFem/FreeFem-sources/tree/develop | defined as a block   |
| /examples/hpddm/stokes-block-2d-PETSc.edp>`__ | system with four     |
|                                               | matrices             |
|                                               | (fieldsplit)         |
+-----------------------------------------------+----------------------+
| `st                                           | Block preconditioner |
| okes-fieldsplit-2d-PETSc.edp <https://github. | (fieldsplit)         |
| com/FreeFem/FreeFem-sources/tree/develop/exam |                      |
| ples/hpddm/stokes-fieldsplit-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `st                                           |                      |
| okes-fieldsplit-3d-PETSc.edp <https://github. |                      |
| com/FreeFem/FreeFem-sources/tree/develop/exam |                      |
| ples/hpddm/stokes-fieldsplit-3d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `maxwell-2d-PETSc.edp <https                  | Direct LU/Cholesky   |
| ://github.com/FreeFem/FreeFem-sources/tree/de |                      |
| velop/examples/hpddm/maxwell-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `maxwell-3d-PETSc.edp <https                  | Multigrid method     |
| ://github.com/FreeFem/FreeFem-sources/tree/de |                      |
| velop/examples/hpddm/maxwell-3d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+

Nonlinear problems
~~~~~~~~~~~~~~~~~~

+-----------------------------------------------+----------------------+
| Filename                                      | Comments             |
|                                               | (preconditioners,    |
|                                               | numerical schemes)   |
+===============================================+======================+
| `bratu-2d-PETSc.edp <htt                      |                      |
| ps://github.com/FreeFem/FreeFem-sources/tree/ |                      |
| develop/examples/hpddm/bratu-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `Newton-2d-PETSc.edp <http                    |                      |
| s://github.com/FreeFem/FreeFem-sources/tree/d |                      |
| evelop/examples/hpddm/Newton-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `                                             | Newton method and    |
| Newton-adaptmesh-2d-PETSc.edp <https://github | *h* adaptivity       |
| .com/FreeFem/FreeFem-sources/tree/develop/exa |                      |
| mples/hpddm/Newton-adaptmesh-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `Newton-vi-2d-PETSc.edp <https:/              | Newton method and a  |
| /github.com/FreeFem/FreeFem-sources/tree/deve | variational          |
| lop/examples/hpddm/Newton-vi-2d-PETSc.edp>`__ | inequality           |
+-----------------------------------------------+----------------------+
| `Newton                                       | Newton method, *h*   |
| -vi-adaptmesh-2d-PETSc.edp <https://github.co | adaptivity, and a    |
| m/FreeFem/FreeFem-sources/tree/develop/exampl | variational          |
| es/hpddm/Newton-vi-adaptmesh-2d-PETSc.edp>`__ | inequality           |
+-----------------------------------------------+----------------------+
| `elasticity-SNES-3d-PETSc.edp <https://githu  | Linear elasiticty    |
| b.com/FreeFem/FreeFem-sources/tree/develop/ex | with a Newton method |
| amples/hpddm/elasticity-SNES-3d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `neo-Hookean-2d-PETSc.edp <https://g          | Nonlinear elasticity |
| ithub.com/FreeFem/FreeFem-sources/tree/develo |                      |
| p/examples/hpddm/neo-Hookean-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `natural-convection-fieldsp                   | Newton method and    |
| lit-2d-PETSc.edp <https://github.com/FreeFem/ | *h* adaptivity       |
| FreeFem-sources/tree/develop/examples/hpddm/n |                      |
| atural-convection-fieldsplit-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `vi-2d-PETSc.edp <                            | Variational          |
| https://github.com/FreeFem/FreeFem-sources/tr | inequalities         |
| ee/develop/examples/hpddm/vi-2d-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+

Time steppers and optimizers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------------------------+----------------------------------+
| Filename                         | Comments (preconditioners,       |
|                                  | numerical schemes)               |
+==================================+==================================+
| `advection-TS-2d-PETSc.edp <http | Implicit and explicit schemes    |
| s://github.com/FreeFem/FreeFem-s |                                  |
| ources/tree/develop/examples/hpd |                                  |
| dm/advection-TS-2d-PETSc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `heat-TS-2d-PETSc.edp            |                                  |
| <https://github.com/FreeFem/Free |                                  |
| Fem-sources/tree/develop/example |                                  |
| s/hpddm/heat-TS-2d-PETSc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `heat-TS-RHS-2d-PETSc.edp <htt   |                                  |
| ps://github.com/FreeFem/FreeFem- |                                  |
| sources/tree/develop/examples/hp |                                  |
| ddm/heat-TS-RHS-2d-PETSc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `minimal-surfa                   | Minimal surface problem          |
| ce-Tao-2d-PETSc.edp <https://git |                                  |
| hub.com/FreeFem/FreeFem-sources/ |                                  |
| tree/develop/examples/hpddm/mini |                                  |
| mal-surface-Tao-2d-PETSc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `orego-Tao-PETSc.edp             |                                  |
|  <https://github.com/FreeFem/Fre |                                  |
| eFem-sources/tree/develop/exampl |                                  |
| es/hpddm/orego-Tao-PETSc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `toy-Tao-PETSc.e                 |                                  |
| dp <https://github.com/FreeFem/F |                                  |
| reeFem-sources/tree/develop/exam |                                  |
| ples/hpddm/toy-Tao-PETSc.edp>`__ |                                  |
+----------------------------------+----------------------------------+

Eigenvalue problems
~~~~~~~~~~~~~~~~~~~

+----------------------------------+----------------------------------+
| Filename                         | Comments (preconditioners,       |
|                                  | numerical schemes)               |
+==================================+==================================+
| `laplace-2d-SLEPc.edp            |                                  |
| <https://github.com/FreeFem/Free |                                  |
| Fem-sources/tree/develop/example |                                  |
| s/hpddm/laplace-2d-SLEPc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `laplace-spherical-harmonics-2   |                                  |
| d-SLEPc.edp <https://github.com/ |                                  |
| FreeFem/FreeFem-sources/tree/dev |                                  |
| elop/examples/hpddm/laplace-sphe |                                  |
| rical-harmonics-2d-SLEPc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `l                               |                                  |
| aplace-torus-2d-SLEPc.edp <https |                                  |
| ://github.com/FreeFem/FreeFem-so |                                  |
| urces/tree/develop/examples/hpdd |                                  |
| m/laplace-torus-2d-SLEPc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `schrodinger-axial-w             |                                  |
| ell-2d-SLEPc.edp <https://github |                                  |
| .com/FreeFem/FreeFem-sources/tre |                                  |
| e/develop/examples/hpddm/schrodi |                                  |
| nger-axial-well-2d-SLEPc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `schro                           |                                  |
| dinger-harmonic-oscillator-1d-SL |                                  |
| EPc.edp <https://github.com/Free |                                  |
| Fem/FreeFem-sources/tree/develop |                                  |
| /examples/hpddm/schrodinger-harm |                                  |
| onic-oscillator-1d-SLEPc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `schro                           |                                  |
| dinger-harmonic-oscillator-2d-SL |                                  |
| EPc.edp <https://github.com/Free |                                  |
| Fem/FreeFem-sources/tree/develop |                                  |
| /examples/hpddm/schrodinger-harm |                                  |
| onic-oscillator-2d-SLEPc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `schrodinger-square-we           |                                  |
| ll-1d-SLEPc.edp <https://github. |                                  |
| com/FreeFem/FreeFem-sources/tree |                                  |
| /develop/examples/hpddm/schrodin |                                  |
| ger-square-well-1d-SLEPc.edp>`__ |                                  |
+----------------------------------+----------------------------------+
| `lapla                           |                                  |
| ce-2d-SLEPc-complex.edp <https:/ |                                  |
| /github.com/FreeFem/FreeFem-sour |                                  |
| ces/tree/develop/examples/hpddm/ |                                  |
| laplace-2d-SLEPc-complex.edp>`__ |                                  |
+----------------------------------+----------------------------------+

Miscellaneous
~~~~~~~~~~~~~

+-----------------------------------------------+----------------------+
| Filename                                      | Comments             |
|                                               | (preconditioners,    |
|                                               | numerical schemes)   |
+===============================================+======================+
| `transpose-solve-PETSc.edp <https://gi        | Solving a transposed |
| thub.com/FreeFem/FreeFem-sources/tree/develop | system               |
| /examples/hpddm/transpose-solve-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `Schur-complement-PETSc.edp <https://git      | Computing an exact   |
| hub.com/FreeFem/FreeFem-sources/tree/develop/ | Schur complement     |
| examples/hpddm/Schur-complement-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `block-PETSc.edp <                            |                      |
| https://github.com/FreeFem/FreeFem-sources/tr |                      |
| ee/develop/examples/hpddm/block-PETSc.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `heat-io-2d.edp                               | Automatic ParaView   |
| <https://github.com/FreeFem/FreeFem-sources/t | animation output     |
| ree/develop/examples/hpddm/heat-io-2d.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `stokes-io-3d.edp <h                          |                      |
| ttps://github.com/FreeFem/FreeFem-sources/tre |                      |
| e/develop/examples/hpddm/stokes-io-3d.edp>`__ |                      |
+-----------------------------------------------+----------------------+
| `buildRecursive.edp <htt                      | Recursive mesh       |
| ps://github.com/FreeFem/FreeFem-sources/tree/ | partitioning (for    |
| develop/examples/hpddm/buildRecursive.edp>`__ | geometric multigrid) |
+-----------------------------------------------+----------------------+
| `withPartitioning.edp <https                  | Connectivity         |
| ://github.com/FreeFem/FreeFem-sources/tree/de | construction with a  |
| velop/examples/hpddm/withPartitioning.edp>`__ | user-supplied        |
|                                               | partitioning         |
+-----------------------------------------------+----------------------+

Reproducible science
~~~~~~~~~~~~~~~~~~~~

+--------------------------------------------------+-------------------+
| Article                                          | Source code       |
+==================================================+===================+
| `Augmented Lagrangian Preconditioner for         | `GitHub           |
| Large-Scale Hydrodynamic Stability               | r                 |
| Analysis <https://www.sciencedire                | epository <https: |
| ct.com/science/article/pii/S0045782519301914>`__ | //github.com/prj- |
|                                                  | /moulin2019al>`__ |
+--------------------------------------------------+-------------------+
| `A Multilevel Schwarz Preconditioner Based on a  | `GitHub           |
| Hierarchy of Robust Coarse                       | repo              |
| Spaces <https://ha                               | sitory <https://g |
| l.archives-ouvertes.fr/hal-02151184/document>`__ | ithub.com/prj-/al |
|                                                  | daas2019multi>`__ |
+--------------------------------------------------+-------------------+
