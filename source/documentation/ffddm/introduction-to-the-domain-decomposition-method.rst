.. _ffddmIntroduction:

Domain Decomposition (DD)
=========================

When the size of a three dimensional problem is large (whatever it means), it is necessary to distribute data among several processors especially for solving linear systems.
A natural way is to do it via domain decomposition.

Mesh Decomposition
------------------

.. raw:: html

   <!--- INSERER FIGURE? The overlap is the minimum width of the overlap between sub-meshes.
   implicit  global renvoie a implicit en fait --->

The starting point is a collection of :math:`N` sub-meshes :math:`(Th_i)_{i=1}^N` that together form a global mesh

.. math:: Th:= \cup_{i=1}^N Th_i\,.

These meshes may be overlapping or not. This decomposition induces a natural decomposition of the global finite element space :math:`Vh` on :math:`Th` into :math:`N` local finite element spaces :math:`(Vh_i)_{i=1}^N` each of them defined on :math:`Th_i`.

**Note** By global, we mean that the corresponding structure can be refered to in the code (most often only) by its local values.
In computer science term, it corresponds to a distributed data where each piece of data is stored by a MPI process.

Distributed Linear Algebra
--------------------------

For a given finite element space :math:`Vh`, the domain decomposition induces a natural decomposition of the set of the global degrees of freedom (d.o.f.) of :math:`Vh` into the :math:`N` subsets of d.o.f.’s :math:`({\mathcal N}_i)_{i=1}^N` each associated with the local finite element space :math:`Vh_i`.
We have thus

.. math:: {\mathcal N} = \cup_{i=1}^N {\mathcal N}_i\,,

but with duplications of some of the d.o.f.’s.

Associated with this decomposition of the set of d.o.f.’s :math:`{\mathcal N}`, a *distributed vector* is a collection of local vectors :math:`({\mathbf V_i})_{1\le i\le N}` so that the values on the duplicated d.o.f.’s are the same.

.. note:: In mathematical terms, it can be described as follows for a real valued problem.
    Let :math:`R_i` be the restriction operator from :math:`\R^{\#{\mathcal N}}` to :math:`\R^{\#{\mathcal N}_i}`, where :math:`\#{\mathcal N}_i` denotes the number of elements of :math:`{\mathcal N}_i`.
    A collection of local vectors :math:`({\mathbf V}_i)_{1\le i\le N}\in \Pi_{i=1}^N \R^{\#{\mathcal N}_i}` is a distributed vector iff there exists a global vector :math:`{\mathbf V}\in\R^{\#{\mathcal N}}` such that for all subset :math:`1\le i\le N`, we have:

    .. math::
        {\mathbf V}_i = R_i\,{\mathbf V}\,.

    We will also say that the collection of local vectors :math:`({\mathbf V}_i)_{1\le i\le N}` is consistent. For a complex valued problem, simply replace :math:`\R` with :math:`\C`.

Partition of Unity Matrices (POUM)
----------------------------------

Let :math:`(D_i)_{1\le i \le N}` be square diagonal matrices of size :math:`\#{\mathcal N}_i` which form a partition of unity in the sense that:

.. math::
     Id_{} = \sum_{i=1}^N R_i^T\,D_i\,R_i\text{ in }\R^{\#{\mathcal N}\times \#{\mathcal N}} \,.

For instance if a degree of freedom is shared by :math:`k` subdomains defining the corresponding entry of the diagonal matrix :math:`D` to be :math:`1/k` yields partition of unity matrices.
The matrices :math:`R_i` and :math:`D_i` are the heart of distributed linear algebra.

Distributed scalar product
~~~~~~~~~~~~~~~~~~~~~~~~~~

For two global vectors :math:`{\mathbf U}` and :math:`{\mathbf V}` of size :math:`\#{\mathcal N}`, the formula for the scalar product :math:`{\mathbf V}^T\,{\mathbf U}=({\mathbf U},\,{\mathbf V})` in terms of their distributed vector counterparts makes use of the partition of unity matrices :math:`(D_i)_{1\le i \le N}` introduced above:

.. math::
   ({\mathbf U}, {\mathbf V}) = \left({\mathbf U}, \sum_{i=1}^N R_i^T D_i R_i {\mathbf V}\right) = \sum_{i=1}^N(R_i {\mathbf U}, D_i R_i {\mathbf V})
   =\sum_{i=1}^N\left({\mathbf U}_i, D_i {\mathbf V}_i\right)\,.

Local scalar products are performed concurrently.
Thus, the implementation is parallel except for the sum which corresponds to a MPI_Reduce call across the :math:`N` MPI processes.
Note also that the implementation relies on the knowledge of a partition of unity so that the FreeFEM syntax is ``dscalprod(Di,u,v)`` or equivalently ``myFEprefix#scalprod(u,v)`` where ``myFEprefix`` is a user defined prefix for the finite element space decomposition, see the :ref:`ffddm documentation <ffddmDocumentationLocalFiniteElementSpaces>`.

.. _ffddmDocumentationUpdate:

Update
~~~~~~

From a collection of local vectors :math:`({\mathbf U}_i)_{1\le i \le N}`, it is possible ensure consistency of the duplicated data by modifying the distributed vector :math:`({\mathbf U}_i)_{1\le i \le N}` by calling the function ``myFEprefix#update(Ui, TRUE)`` where ``myFEprefix`` is the user defined prefix that refers to the finite element space decomposition.
This function performs the following operation for all :math:`1\le i \le N`:

.. math::
    {\mathbf U}_i \leftarrow R_i\, \sum_{j=1}^N R_j^T D_j {\mathbf U}_j

.. note:: The implementation corresponds to

    .. math::
        {\mathbf U}_i \leftarrow R_i \sum_{j=1}^N R_j^T D_j {\mathbf U}_j = D_i {\mathbf U}_i + \sum_{j\in \mathcal{O}(i)} R_i\,R_j^T\,D_j {\mathbf U}_j

    where :math:`\mathcal{O}(i)` is the set of neighbors of subdomain :math:`i`.
    Therefore, the matrix vector product is computed in three steps: 
    
    - concurrent computing of :math:`D_j {\mathbf U}_j` for all :math:`1\le j\le N`; 
    - neighbor to neighbor MPI-communications from subdomain :math:`j` to subdomain :math:`i`  (:math:`R_i\,R_j^T`) ; 
    - concurrent sum of neighbor contributions.

Distributed Matrix and Vector resulting from a variational formulation
----------------------------------------------------------------------

The discretization of a variational formulation on the global mesh :math:`Th` yields a global matrix :math:`A` and a global right hand side :math:`\mathbf{RHS}`.
Thanks to the sparsity of finite element matrices for partial differential equations and thanks to the overlap between subdomains, the knowledge of the local matrix :math:`R_i A R_i^T` on each subdomain :math:`1\le i\le N` is sufficient to perform the matrix-vector product :math:`A\times \mathbf{U}` for any global vector :math:`\mathbf{U}`.
Once the problem has been set up by a call to ``ffddmsetupOperator(myprefix, myFEprefix, myVarf)``, the matrix-vector product is performed by calling the function ``myprefix#A(Ui)`` where ``myprefix`` is a user defined prefix that refers to the problem at hand which itself implicitly refers to the triplet (domain decomposition, finite element, variational formulation).
See more on problem definition in this :ref:`documentation <ffddmDocumentationDefineProblemToSolve>` and more on distributed linear algebra in chapter 8 of `"An Introduction to Domain Decomposition Methods: algorithms, theory and parallel implementation" SIAM 2015 <http://bookstore.siam.org/ot144/>`__.

Distributed Linear Solvers
--------------------------

In many cases, we are interested in the solution of the problem in terms of the vector of d.o.f.’s :math:`\mathbf{X}` that satisfies:

.. math:: A\, \mathbf{X} = \mathbf{RHS}\,.

``ffddm`` offers two parallel solvers: :ref:`direct factorization <ffddmIntroductionDisitributedDirectSolvers>` and iterative preconditioned solvers via :ref:`Schwarz <ffddmIntroductionSchwarzMethods>` domain decomposition methods.

.. _ffddmIntroductionDisitributedDirectSolvers:

Distributed Direct Solvers
~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to benefit from the sparsity of the matrix arising from a finite element discretization of a partial differential equation, a variant of Gauss elimination, the frontal method, that automatically avoids a large number of operations involving zero terms was developed.
A frontal solver builds a :math:`LU` or Cholesky decomposition of a sparse matrix given as the assembly of element matrices by eliminating equations only on a subset of elements at a time.
This subset is called the *front* and it is essentially the transition region between the part of the system already finished and the part not touched yet.
These methods are basically sequential since the unknowns are processed the one after another or one front after another.
In order to benefit from multicore processors, a `multifrontal solver <https://en.wikipedia.org/wiki/Multifrontal_method>`__ is an improvement of the frontal solver that uses several independent fronts at the same time.
The fronts can be worked on by different processors, which enables parallel computing. ``ffddm`` provides an interface to the parallel sparse direct solver `MUMPS <http://mumps.enseeiht.fr/>`__. These methods have the advantage to be very robust and to have a predictable cost. The main drawback is the memory requirement which can be prohibitive especially for three-dimensional problems. 

.. _ffddmIntroductionSchwarzMethods:

Schwarz methods
~~~~~~~~~~~~~~~

These methods are part of the large family of preconditioned iterative solvers. When considering the solve of the equation :math:`A\, \mathbf{X} = \mathbf{RHS}`, a preconditioner is a linear operator that approximates the inverse of :math:`A` and whose cost of the associated matrix-vector product is much cheaper than solving the original linear system. It enables to accelerate the solution of the latter with Krylov type methods such as the conjugate gradient (in the symmetric positive definite case), GMRES or BiCGSTAB in the general case. Two options are possible. 

Left preconditioning: the preconditioner is applied to the left of the equation 

.. math::
   M^{-1}  A\, \mathbf{X} =  M^{-1} \mathbf{RHS}\,.

and the Krylov method is applied to the left preconditioned system with a residual that is preconditioner dependent. 

Right preconditioning: the preconditioner is inserted on the right of the operator:

.. math::
    A\, M^{-1}  \mathbf{Y} =  \mathbf{RHS}\, \text{ where } \mathbf{X} =  M^{-1}  \mathbf{Y}.

and the Krylov method is applied to the right preconditioned system with a residual that is preconditioner independent.  

In both cases, if the preconditioner is efficient the number of iterations to get a converged solution is much smaller than the number of iterations of the Krylov method applied to the original equation :math:`A\, \mathbf{X} = \mathbf{RHS}`.  Although right preconditioning seems more intricate, it is much safer to use since the convergence is checked on a residual that does not depend on the preconditioner.

In the sequel, we consider the solution of the equation :math:`A\, \mathbf{X} = \mathbf{RHS}` preconditioned by domain decomposition methods and with a **flexible GMRES** Krylov method which is thus necessarily right preconditioned. 

Restricted Additive Schwarz (RAS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The RAS preconditioner reads:

.. math::
   M^{-1}_{RAS} := \sum_{j=1}^N R_j^T D_j (R_j\, A\,R_j^T)^{-1} R_j\,,

where for each subdomain :math:`j` the restriction matrix :math:`R_j` and  the partition of unity matrix :math:`D_j`have been introduced above. Note that in the original ASM (additive Schwarz method) preconditioner the partition of unity is dropped. The application of the operator :math:`M^{-1}_{RAS}` to a global right hand side :math:`\mathbf{RHS}` is detailed below. Recall that this global vector is distributed among processes via the local vectors :math:`(\mathbf{RHS}_i)_{i=1}^N`. Let :math:`A_{j}` denote the local matrix :math:`(R_j\, A\,R_j^T)`. The local vector in subdomain :math:`i` resulting from the matrix vector product :math:`M^{-1}_{RAS}\, \mathbf{RHS}` consists in computing:

.. math::
   R_i\, \sum_{j=1}^N R_j^T\,D_j\, A_{j}^{-1}\,\, \mathbf{ RHS}_j
   = D_i\, A_{i}^{-1}\, \mathbf{ RHS}_i + \sum_{j\in \mathcal{O}(i)} (R_i\,R_j^T)\,D_j\, A_{j}^{-1}\, \mathbf{ RHS}_j\,.

This task is performed by first solving concurrently on all subdomains a linear system for :math:`{\mathbf Y}_j` for all :math:`1\le j \le N`:

.. math::
   A_{j}\, {\mathbf Y}_j = \mathbf{RHS}_j\,.

Each local vector :math:`{\mathbf Y}_j` is weighted by the partition of unity matrix :math:`D_j`.
Then data transfers between neighboring subdomains implement the :math:`R_i\,R_j^T\,D_j\,{\mathbf Y}_j` formula.
The contribution from neighboring subdomains are summed locally. This
pattern is very similar to that of the :ref:`update <ffddmDocumentationUpdate>` procedure.

Optimized Restricted Additive Schwarz (ORAS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ORAS preconditioner may be seen as a variant of the RAS preconditioner.
It reads:

.. math::
   M^{-1}_{RAS} := \sum_{j=1}^N R_j^T D_j\, B_j^{-1}\, R_j\,

where :math:`B_j` are local matrices of size :math:`\#{\mathcal N}_j \times \#{\mathcal N}_j` for :math:`1\le j \le N`.
This variant is very useful when dealing with wave propagation phenomena such as Helmholtz problems in acoustics or Maxwell system in the frequency domain for electromagnetism.
Defining :math:`B_j` as the discretization of the physical equation with impedance conditions on the boundary of the subdomain has been proved to be a good choice.

Two level methods
^^^^^^^^^^^^^^^^^

The RAS and ORAS methods are called a one-level method in the sense that sub-domains only interact with their direct neighbors. For some problems such as Darcy problems or static elasticity problems and when the number of subdomains is large, such one-level methods may suffer from a slow convergence.
The fix is to add to the preconditioner an auxiliary coarse problem that couples all subdomains at each iteration and is inexpensive to calculate.

In mathematical terms, we first choose  a full rank rectangular matrix  :math:`Z\in\R^{\#{\mathcal N}\times NC}` where :math:`NC \ll \#{\mathcal N}` denotes the dimension of the coarse space spanned by the columns of :math:`Z`. We also pick a coarse matrix :math:`A_C\in \R^{N_C\times N_C}`. A generic one-level method preconditioner :math:`M_1^{-1}` is enriched by a solve on the coarse space. The simplest correction formula is additive:

.. math::
  M_2^{-1} := Z \,A_C^{-1}\,Z^T + M_1^{-1}

Other correction formulas are given in :ref:`documentation <ffddmDocumentationTwoLevelPreconditioners>`.

We consider two ways to build :math:`Z` and thus the coarse space and the coarse problem :math:`A_C`, see below :ref:`Coarse Mesh <ffddmIntroductionCoarseMesh>` and :ref:`GenEO <ffddmIntroductionGeneo>`

.. _ffddmIntroductionCoarseMesh:

Coarse Mesh
'''''''''''

A first possibility is to discretize the problem on a coarse mesh, following the same principle as multi-grid methods.
For 3-D problems, a coarsening of the mesh size by a factor 2, reduces by a factor :math:`2^3=8` the size of the coarse problem which is then easier to solve by a direct method. Then, :math:`Z` is the interpolation matrix from the coarse finite element space to the fine one.


.. _ffddmIntroductionGeneo:

GenEO
'''''

For highly heterogeneous or anisotropic problems, two level methods based on coarse meshes might fail and a more sophisticated construction must be used.
A provable robust coarse space called GenEO is built by first solving the following local generalized eigenvalue problem in parallel for each subdomain :math:`1\le i\le N`, where :math:`A_i^{\text{Neu}}` denotes the local matrix resulting from the variational formulation:

.. math::
   D_i A_i D_i\, V_{i,k} = \lambda_{i,k}\, A_i^{\text{Neu}} \,V_{i,k}

The eigenvectors selected to enter the coarse space correspond to eigenvalues :math:`\lambda_{i,k} \ge \tau`, where the threshold parameter :math:`\tau` is user-defined.
The precise formulas are given in this :ref:`documentation <ffddmDocumentationBuildingGeneoCoarseSpace>`.
From a mathematical point of view, it has been proved that for a symmetric positive definite matrix :math:`A`, the spectrum of the preconditioned by the two-level method with a GenEO coarse space lies in the interval :math:`[\displaystyle \frac{1}{1+k_1\,\tau} , k_0 ]`.

**Note** A heuristic that justifies this construction is as follows.
We first introduce the Additive Schwarz method (ASM) which can be seen as a symmetrized variant of the RAS preconditioner:

.. math::
       M_{ASM}^{-1} := \sum_{j=1}^N R_j^T A_j^{-1} R_j\,.

It can be proved that the lower bound for the eigenvalue of :math:`M_{ASM}^{-1}\,A` is close to zero (which is bad for convergence) whereas the upper bound depends only on the number of neigbors of a subdomain (which is good for convergence).

Second, we also introduce the following preconditioner :math:`M^{-1}_{NN}`:

.. math::
       M^{-1}_{NN} := \sum_{1\le j\le N} D_i\,(A_j^{\text{Neu}})^{-1} D_j\,.

We have a very good lower bound for the preconditioned operator :math:`M^{-1}_{NN}\,A` that does not depend on the number of subdomains but only on the maximum multiplicity of intersections :math:`k_1` (which is good for convergence).
But the upper bound for this preconditioner is very large (which is bad for convergence).

Now, if we compare formulas for :math:`M^{-1}_{NN}` and :math:`M^{-1}_{ASM}`, we may suspect that vectors :math:`\mathbf{V}_{ik}` for which :math:`D_i\, (A_i^{\text{Neu}})^{-1}\,D_i\,\mathbf{V}_{ik}` and :math:`A_{i}^{-1}\,\mathbf{V}_{ik}` have very different values are responsible for the slow convergence and should contribute to the coarse space.
This is a way to interpret the above generalized eigenvalue problem which controls the lower bound of the two-level preconditioned system.
