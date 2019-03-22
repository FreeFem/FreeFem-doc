.. _ffddmExamples:

Examples
========

+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| File name                             | :math:`M^{-1}_1` | :math:`M^{-1}_2` | inexact CS  | comments                                 |
+=======================================+==================+==================+=============+==========================================+
| `diffusion-3d-minimal-direct.edp`_    |                  |                  |             | direct solver MUMPS                      |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `diffusion-3d-minimal-ddm.edp`_       | RAS              | GenEO            |             |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `diffusion-3d-simple.edp`_            | RAS              | GenEO            |             | comparison with direct solver            |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `diffusion-2d-thirdlevelgeneo.edp`_   | RAS              | GenEO            | RAS + GenEO |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `elasticity-3d-simple.edp`_           | RAS              | GenEO            |             |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `elasticity-3d-thirdlevelgeneo.edp`_  | RAS              | GenEO            | RAS + GenEO |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `Helmholtz-2d-simple.edp`_            | ORAS             | Coarse Mesh      |             |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `Helmholtz-2d-marmousi.edp`_          | ORAS             | Coarse Mesh      |             |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `Helmholtz-3d-simple.edp`_            | ORAS             | Coarse Mesh      |             |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `Helmholtz-3d-overthrust.edp`_        | ORAS             |                  |             |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `Helmholtz-2d-HPDDM-BGMRES.edp`_      | ORAS             |                  |             | multi-rhs Block GMRES with HPDDM         |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `Navier-2d-marmousi2.edp`_            | ORAS             | Coarse Mesh      |             |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `Maxwell-3d-simple.edp`_              | ORAS             | Coarse Mesh      |             |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `Maxwell_Cobracavity.edp`_            | ORAS             | Coarse Mesh      | ORAS        |                                          |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `natural_convection.edp`_             | ORAS             | Coarse Mesh      |             | nonlinear                                |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `natural_convection_3D_obstacle.edp`_ | ORAS             | Coarse Mesh      |             | nonlinear                                |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+
| `Richards-2d.edp`_                    | RAS              |                  |             | nonlinear time dependent mesh adaptation |
+---------------------------------------+------------------+------------------+-------------+------------------------------------------+

.. _diffusion-3d-minimal-direct.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/diffusion-3d-minimal-direct.edp
.. _diffusion-3d-minimal-ddm.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/diffusion-3d-minimal-ddm.edp
.. _diffusion-3d-simple.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/diffusion-3d-simple.edp
.. _diffusion-2d-thirdlevelgeneo.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/diffusion-2d-thirdlevelgeneo.edp
.. _elasticity-3d-simple.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/elasticity-3d-simple.edp
.. _elasticity-3d-thirdlevelgeneo.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/elasticity-3d-thirdlevelgeneo.edp
.. _Helmholtz-2d-simple.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-2d-simple.edp
.. _Helmholtz-2d-marmousi.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-2d-marmousi.edp
.. _Helmholtz-3d-simple.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-3d-simple.edp
.. _Helmholtz-3d-overthrust.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-3d-overthrust.edp
.. _Helmholtz-2d-HPDDM-BGMRES.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-2d-HPDDM-BGMRES.edp
.. _Navier-2d-marmousi2.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Navier-2d-marmousi2.edp
.. _Maxwell-3d-simple.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Maxwell-3d-simple.edp
.. _Maxwell_Cobracavity.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Maxwell_Cobracavity.edp
.. _natural_convection.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/natural_convection.edp
.. _natural_convection_3D_obstacle.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/natural_convection_3D_obstacle.edp
.. _Richards-2d.edp: https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Richards-2d.edp
