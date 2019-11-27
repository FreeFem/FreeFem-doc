PETSc and SLEPc
===============

FreeFEM is interfaced with `PETSc <https://www.mcs.anl.gov/petsc/>`__ and `SLEPc <http://slepc.upv.es/>`__ which offer a wide range of sequential or parallel linear or nonlinear solvers, time steppers, optimizers, and eigensolvers. In particular, it gives access transparently (without much changes to user code) to: distributed and multithreaded direct solvers (`PARDISO <https://software.intel.com/en-us/mkl-developer-reference-fortran-intel-mkl-pardiso-parallel-direct-sparse-solver-interface>`__, `MUMPS <http://mumps.enseeiht.fr/>`__, `SuperLU <https://crd.lbl.gov/departments/applied-mathematics/scalable-solvers/research-areas/superlu-dist/>`__), multigrid solvers (`hypre <https://hypre.readthedocs.io/en/latest/>`__, GAMG), domain decomposition methods (block Jacobi, ASM, `HPDDM <https://github.com/hpddm/hpddm/>`__). For a detailed introduction to these tools, interested readers are referred to the tutorial `Introduction to FreeFEM with an emphasis on parallel computing <http://jolivet.perso.enseeiht.fr/FreeFem-tutorial/>`__.

In most of the scripts listed below, the following standard procedure is used.
        * Load an initial sequential mesh (in 2D or 3D).

        * Partition the mesh and generate connectivity information according to the number of processes.

        * Provide these information to PETSc so that subsequent computations may be done in a distributed fashion.

Combining the power and flexibility of PETSc with the ease-of-use of FreeFEM may help design multiphysics solvers, e.g., for `Navier--Stokes equations <https://www.sciencedirect.com/science/article/pii/S0045782519301914>`__, advanced matrix-free discretizations, and such.

.. toctree::

  examples
