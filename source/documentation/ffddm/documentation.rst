.. role:: freefem(code)
   :language: freefem

.. _ffddmDocumentation:

ffddm documentation
===================

Minimal example
---------------

.. code-block:: freefem
   :linenos:

   macro dimension 3// EOM            // 2D or 3D

   include "ffddm.idp"

   int[int] LL = [2,2, 1,2, 2,2];
   mesh3 ThGlobal = cube(10, 10, 10, [x, y, z], label = LL);      // global mesh

   macro grad(u) [dx(u), dy(u), dz(u)]// EOM    // three-dimensional gradient

   macro Varf(varfName, meshName, VhName)
       varf varfName(u,v) = int3d(meshName)(grad(u)'* grad(v)) + int3d(meshName)(v) + on(1, u = 1.0);
   // EOM

   // Domain decomposition
   ffddmbuildDmesh( Lap , ThGlobal , mpiCommWorld )

   macro def(i)i// EOM                         // scalar field definition
   macro init(i)i// EOM                        // scalar field initialization
   ffddmbuildDfespace( Lap , Lap , real , def , init , P1 )

   ffddmsetupOperator( Lap ,Lap , Varf )

   real[int] rhsi(0);
   ffddmbuildrhs( Lap , Varf , rhsi )

   LapVhi def(ui);

   //Direct solve
   ui[] = Lapdirectsolve(rhsi);

   Lapwritesummary

   ffddmplot(Lap,ui,"u");

.. _ffddmDocumentationOverlappingMeshDecomposition:

Overlapping mesh decomposition
------------------------------

.. code-block:: freefem
   :linenos:

   ffddmbuildDmesh(pr,Th,comm)

decomposes the mesh **Th** into overlapping submeshes.
The mesh will be distributed over the mpi ranks of communicator **comm**.
This will create and expose variables whose names will be prefixed by **pr**, see below (# is the concatenation operator).
The way the initial mesh **Th** is partitioned depends on the value of :ref:`ffddmpartitioner <ffddmParametersGlobal>`.

The size of the overlap between subdomains (its width in terms of number of mesh elements) is given by :ref:`ffddmoverlap <ffddmParametersGlobal>`.

The level of refinement of the resulting submeshes with respect to the input mesh **Th** is given by :ref:`ffddmsplit <ffddmParametersGlobal>`.

If :ref:`ffddmexclude <ffddmParametersGlobal>` :math:`\neq 0`, the first :ref:`ffddmpCS <ffddmParametersGlobal>` mpi ranks of **comm** will be excluded from the spatial domain decomposition, in order to dedicate them later to the coarse problem (for two-level preconditioners).

The label of the new border of the submeshes (the interface between the subdomains) is given by :ref:`ffddminterfacelabel <ffddmParametersGlobal>`.

**defines**:

-  ``int pr#npart`` number of subdomains for this decomposition; should be equal to mpiSize(\ **comm**) - :ref:`ffddmexclude <ffddmParametersGlobal>` * :ref:`ffddmpCS <ffddmParametersGlobal>`
-  ``meshN[int] pr#aTh`` array (size ``pr#npart``) of local meshes of the subdomains.
   In the standard parallel case, only the local mesh for this mpi rank ``pr#aTh[mpiRank(pr#commddm)]`` is defined (unless this mpi rank is excluded from the spatial domain decomposition, i.e. ``prmesh#excluded`` = 1, see below).
   In the sequential case, all local meshes are defined.
-  ``meshN pr#Thi`` the local mesh of the subdomain for this mpi rank, i. e. ``pr#aTh[mpiRank(pr#commddm)]`` in the parallel case -  ``int pr#numberIntersection`` the number of neighbors for this mpi rank
-  ``int[int] pr#arrayIntersection`` the list of neighbor ranks in ``pr#commddm`` for this mpi rank
-  ``int pr#pCS`` equal to :ref:`ffddmpCS <ffddmParametersGlobal>`
-  ``int pr#exclude`` equal to :ref:`ffddmexclude <ffddmParametersGlobal>`
-  ``int pr#excluded`` *true* if :ref:`ffddmexclude <ffddmParametersGlobal>` is *true* (:math:`\neq 0`) and mpiRank(\ **comm**) < ``pr#pCS``.
   In this case, this mpi rank will be excluded from the spatial domain decomposition and will only work on the coarse problem.
-  ``mpiComm pr#commddm`` mpi communicator for ranks participating in the spatial domain decomposition (ranks 0 to ``pr#npart``-1 in **comm** if ``pr#exclude`` is *false*, ranks ``pr#pCS`` to ``pr#pCS``\ +\ ``pr#npart``-1 otherwise)
-  ``mpiComm pr#commCS`` mpi communicator for ranks participating in the assembly and resolution of the coarse problem for two-level preconditioners (ranks 0 to ``pr#pCS`` - 1 in **comm**)
-  ``mpiComm pr#commself`` self mpi communicator (this mpi rank only), used for factorizing local matrices

.. raw:: html

   <!--
   ***For advanced users***:

   - `int pr#binexactCS`
   - `int pr#inexactCSsplit`
   - `int pr#isincomm`
   - `meshN[int] pr#aThborder`
   -->

**Remark for sequential use** (see :ref:`-seqddm <ffddmParametersCommandLine>`):
    - ``meshN[int] pr#aTh`` array (size ``pr#npart``) of local meshes of the subdomains

.. raw:: html

   <!--
    int pr#binexactgeneoCS

   fespace pr#VhiP1(pr#Thi,P1);

   pr#VhiP1[int] pr#partitionIntersectionbasei(0);

   meshN pr#Thglob = minimalMesh;

   matrix[int] pr#RihP1(pr#npart);
   pr#VhiP1[int] pr#DP1(pr#npart);

   NewMacro pr#mpicomm()comm EndMacro

   ***depends on***:
   - [ffddmpartitioner](parameters.md#global-parameters)
   - [ffddmpCS](parameters.md#global-parameters)
   - [ffddmexclude](parameters.md#global-parameters)
   - [ffddmoverlap](parameters.md#global-parameters)
   - [ffddmsplit](parameters.md#global-parameters)
   - [ffddminterfacelabel](parameters.md#global-parameters)

   ***see also***:
   -->

.. _ffddmDocumentationLocalFiniteElementSpaces:

Local finite element spaces
---------------------------

.. code-block:: freefem
   :linenos:

   ffddmbuildDfespace(pr,prmesh,scalar,def,init,Pk)

builds the local finite element spaces and associated distributed operators on top of the mesh decomposition **prmesh**.
This will create and expose variables whose names will be prefixed by **pr**, see below.
It is assumed that :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>` has already been called with prefix **prmesh** in order to build the mesh decomposition.

The local finite element spaces of type **Pk** (where **Pk** is the type of finite element: P1, [P2,P2,P1], …) are defined on the local meshes of the subdomains based on the mesh decomposition previously created with prefix **prmesh**.

**scalar** determines the type of data for this finite element: *real* or *complex*.

Two macros, **def** and **init**, are needed: **def** specifies how to define a finite element function in the finite element space **Pk**, and **init** specifies how to interpolate a scalar function onto the (possibly multiple) components of **Pk**. Two examples are given below:

For scalar P2 finite elements and complex-valued problems:

.. code-block:: freefem
   :linenos:

   macro def(u) u// EOM
   macro init(u) u// EOM
   ffddmbuildDfespace(myFEprefix,mymeshprefix,complex,def,init,P2)

For vectorial [P2,P2,P1] finite elements and real-valued problems:

.. code-block:: freefem
   :linenos:

   macro def(u) [u, u#B, u#C]// EOM
   macro init(u) [u, u, u]// EOM
   ffddmbuildDfespace(myFEprefix,mymeshprefix,real,def,init,[P2,P2,P1])

In practice, this builds the necessary distributed operators associated to the finite element space: the local partition of unity functions :math:`(D_i)_{i=1,...,N}` (see ``pr#Dk`` and ``pr#Dih`` below) as well as the function ``pr#update`` (see below) which synchronizes local vectors :math:`(u_i)_{i=1,...,N}` between neighboring subdomains, performing the equivalent of :math:`u_i = R_i (\sum_{j=1}^N R_j^T u_j)` or :math:`u_i = R_i (\sum_{j=1}^N R_j^T D_j u_j)` in a distributed parallel environment.

``pr#scalprod`` (see below) performs the parallel scalar product for vectors defined on this finite element.

**defines**:

-  ``pr#prmesh`` macro, saves the parent prefix **prmesh** of the mesh decomposition
-  ``pr#K`` macro, saves the type of data **scalar** for this finite element space (*real* or *complex*)
-  ``func pr#fPk`` saves the type of finite element **Pk**, e.g. \ *P1, [P2,P2,P1], …*
-  ``fespace pr#Vhi`` the local finite element space for this mpi rank, defined on the local mesh ``prmesh#Thi``
-  ``int pr#Ndofglobal`` the total number of degrees of freedom :math:`n` for this finite element discretization
-  ``pr#mdef`` macro, saves the macro **def** giving the definition of a finite element function in the finite element space **Pk**
-  ``pr#minit`` macro, saves the macro **init** specifying how to interpolate a scalar function onto the (possibly multiple) components of a finite element function of **Pk**.
   This is used to create the local partition of unity function in ``pr#Vhi``, by interpolating the local P1 partition of unity function onto the components of ``pr#Vhi``.
   For non Lagrange finite element spaces (e.g. *RT0*, *Edge03d*, …), see :ref:`ffddmbuildDfespaceEdge <ffddmDocumentationPartitionUnityEdge>`.
-  ``pr#K[int][int] pr#Dk`` array (size ``prmesh#npart``) of local partition of unity vectors in the subdomains, equivalent to :math:`(D_i)_{i=1,...,N}`.
   In the standard parallel case, only the local partition of unity vector for this mpi rank ``pr#Dk[mpiRank(prmesh#commddm)]`` is defined (unless this mpi rank is excluded from the spatial domain decomposition, i. e. ``prmesh#excluded`` = 1).
   In the sequential case, all local partition of unity vectors are defined.
-  ``matrix<pr#K>[int] pr#Dih`` array (size ``prmesh#npart``) similar to ``pr#Dk`` but in *matrix* form, allowing for easier *matrix*-*matrix* multiplications.
   ``pr#Dih[i]`` is a diagonal matrix, with the diagonal equal to ``pr#Dk[i]``.
-  ``fespace pr#Vhglob`` the global finite element space defined on the global mesh ``prmesh#Thglob``.
   Defined only if :ref:`-noGlob <ffddmParametersCommandLine>` is not used.
-  ``matrix<pr#K>[int] pr#Rih`` array (size ``prmesh#npart``) of restriction matrices from the global finite element space to the local finite element spaces on the local submeshes of the subdomains.
   In the standard parallel case, only the restriction matrix for this mpi rank ``pr#Rih[mpiRank(prmesh#commddm)]`` is defined (unless this mpi rank is excluded from the spatial domain decomposition, i. e. ``prmesh#excluded`` = 1).
   In the sequential case, all restriction matrices are defined. The restriction matrices ``pr#Rih`` are defined only if :ref:`-noGlob <ffddmParametersCommandLine>` is not used.
-  ``func int pr#update(scalar[int] ui, bool scale)`` The function ``pr#update`` synchronizes the local vector *ui* between subdomains by exchanging the values of *ui* shared with neighboring subdomains (in the overlap region) using point-to-point MPI communications.
   If *scale* is *true*, *ui* is multiplied by the local partition of unity beforehand.
   This is equivalent to :math:`u_i = R_i (\sum_{j=1}^N R_j^T u_j)` when *scale* is *false* and :math:`u_i = R_i (\sum_{j=1}^N R_j^T D_j u_j)` when *scale* is *true*.
-  ``func scalar pr#scalprod(scalar[int] ai, scalar[int] bi)`` The function ``pr#scalprod`` computes the global scalar product of two vectors whose local restriction to the subdomain of this mpi rank are *ai* and *bi*.
   The result is computed as :math:`\sum_{j=1}^N (D_j a_j, b_j)`.

.. raw:: html

   <!--
   ***Remark:***


   ***For advanced users***:

   matrix<pr#K>[int] pr#restrictionIntersection(0);

   NewMacro pr#mdefpart udefpart EndMacro

   NewMacro pr#minitpart uinitpart EndMacro

   func pr#fPkP0 = mPkP0;

   pr#K[int][int] pr#rcv(0);
   pr#K[int][int] pr#snd(0);

   ***depends on***:

   ***see also***:

   - **[`ffddmbuildDfespaceEdge`](#local-finite-element-spaces)**
   -->

.. _ffddmDocumentationDefineProblemToSolve:

Define the problem to solve
---------------------------

.. code-block:: freefem
   :linenos:

   ffddmsetupOperator(pr,prfe,Varf)

builds the distributed operator associated to the variational problem given by **Varf**, on top of the distributed finite element space **prfe**.
This will create and expose variables whose names will be prefixed by **pr**, see below.
It is assumed that :ref:`ffddmbuildDfespace <ffddmDocumentationLocalFiniteElementSpaces>` has already been called with prefix **prfe** in order to define the distributed finite element space.

In practice, this builds the so-called local ‘Dirichlet’ matrices :math:`A_i = R_i A R_i^T`, the restrictions of the global operator :math:`A` to the subdomains (see ``pr#aRd``\ below).
The matrices correspond to the discretization of the bilinear form given by the macro **Varf**, which represents the abstract variational form of the problem.
These matrices are then used to implement the action of the global operator :math:`A` on a local vector (the parallel matrix-vector product with :math:`A`), see ``pr#A`` below.

At this point, we already have the necessary data to be able to solve the problem with a parallel direct solver (*MUMPS*), which is the purpose of the function ``pr#directsolve`` (see below).
See :ref:`ffddmbuildrhs <ffddmDocumentationBuildRhs>` for building the right-hand side.

The macro **Varf** is required to have three parameters: the name of the variational form, the mesh, and the finite element space.
The variational form given in this ‘abstract’ format will then be used by *ffddm* to assemble the discrete operators by setting the appropriate mesh and finite element space as parameters.
An example is given below:

.. code-block:: freefem
   :linenos:

   macro myVarf(varfName, meshName, VhName)
       varf varfName(u,v) = int3d(meshName)(grad(u)''* grad(v)) + on(1, u = 1.0);
   // EOM

   ffddmsetupOperator(myprefix,myFEprefix,myVarf)

**Remark** In this simple example, the third parameter *VhName* is not used.
However, for more complex cases such as non-linear or time dependent problems where the problem depends on a solution computed at a previous step, it is useful to know for which discrete finite element space the variational form is being used.
See for example TODO

**defines**:

-  ``pr#prfe`` macro, saves the parent prefix **prfe** of the finite element space
-  ``int pr#verbosity`` the level of verbosity for this problem, initialized with the value of :ref:`ffddmverbosity <ffddmParametersGlobal>`
-  ``pr#writesummary`` macro, prints a summary of timings for this problem, such as the time spent to assemble local matrices or solve the linear system.
-  ``matrix<prfe#K> pr#Aglobal`` the global matrix :math:`A` corresponding to the discretization of the variational form given by the macro **Varf** on the global finite element space ``prfe#Vhglob``.
   Defined only in the sequential case.
-  ``matrix<prfe#K>[int] pr#aRd`` array (size ``prfe#prmesh#npart``) of so-called local ‘Dirichlet’ matrices in the subdomains; these are the restrictions of the global operator to the subdomains, equivalent to :math:`A_i = R_i A R_i^T` with :math:`A` the global matrix corresponding to the discretization of the variational form given by the macro **Varf** on the global finite element space.
   In the standard parallel case, only the local matrix for this mpi rank ``pr#aRd[mpiRank(prmesh#commddm)]`` is defined (unless this mpi rank is excluded from the spatial domain decomposition, i. e. ``prmesh#excluded`` = 1).
   In the sequential case, all local matrices are defined.
-  ``func prfe#K[int] pr#A(prfe#K[int] &ui)`` The function ``pr#A`` computes the parallel matrix-vector product, i.e. the action of the global operator :math:`A` on the local vector :math:`u_i`.
   The computation is equivalent to :math:`R_i (\sum_{j=1}^N R_j^T D_j A_j u_j)` and is performed in parallel using local matrices ``pr#aRd`` and the function ``prfe#update``.
   In the sequential case, the global matrix ``pr#Aglobal`` is used instead.
-  ``func prfe#K[int] pr#AT(prfe#K[int] &ui)`` Similarly to ``pr#A``, The function ``pr#AT`` computes the action of :math:`A^T`, the transpose of the global operator :math:`A`, on `u_i`.
-  ``func prfe#K[int] pr#directsolve(prfe#K[int]& rhsi)`` The function ``pr#directsolve`` allows to solve the linear system :math:`A x = b` in parallel using the parallel direct solver *MUMPS*.
   The matrix is given to *MUMPS* in distributed form through the local matrices ``pr#aRd``.
   The input *rhsi* is given as a distributed vector (*rhsi* is the restriction of the global right-hand side :math:`b` to the subdomain of this mpi rank, see :ref:`ffddmbuildrhs <ffddmDocumentationBuildRhs>`) and the returned vector is local as well.

**Remark: rectangular operators**

It is possible to define a non-square distributed operator where the variational form takes two different finite element spaces of unknown and test functions. This is done through macro **ffddmsetupOperatorRect** which takes two FE prefixes (which must be defined on the same mesh prefix), see below:

.. code-block:: freefem
   :linenos:

   macro myVarf(varfName, meshName, VhName)
       varf varfName([u, uB, uC], [q]) = int3d(meshName)(div(u) * q);
   // EOM

   ffddmsetupOperatorRect(myprefix,myFEprefixV,myFEprefixP,myVarf)

.. raw:: html

   <!--
   NewMacro pr#plot(u,s)

   ***For advanced users***:

   NewMacro pr#fromVhi(ui,VhName,res)

   ***depends on***:

   - [ffddmverbosity](parameters.md#global-parameters)
   -->

--------------

.. _ffddmDocumentationBuildRhs:

.. code-block:: freefem
   :linenos:

   ffddmbuildrhs(pr,Varfrhs,rhs)

builds the right-hand side associated to the variational form given by **Varfrhs** for the problem corresponding to prefix **pr**.
The resulting right-hand side vector **rhs** corresponds to the discretization of the abstract linear form given by the macro **Varfrhs** (see :ref:`ffddmsetupOperator <ffddmDocumentationDefineProblemToSolve>` for more details on how to define the abstract variational form as a macro).

The input vector **rhs** is resized and contains the resulting local right-hand side :math:`R_i b`, the restriction of the global right-hand side :math:`b` to the subdomain of this mpi rank.
In the sequential case, the global right-hand side vector :math:`b` is assembled instead.

An example is given below:

.. code-block:: freefem
   :linenos:

   macro myVarfrhs(varfName, meshName, VhName)
       varf varfName(u,v) = intN(meshName)(v) + on(1, u = 1.0);
   // EOM

   real[int] rhsi(0);
   ffddmbuildrhs(myprefix,myVarfrhs,rhsi)

.. _ffddmDocumentationOneLevelPreconditioners:

One level preconditioners
-------------------------

.. code-block:: freefem
   :linenos:

   ffddmsetupPrecond(pr,VarfPrec)

builds the one level preconditioner for problem **pr**.
This will create and expose variables whose names will be prefixed by **pr**, see below.
It is assumed that :ref:`ffddmsetupOperator <ffddmDocumentationDefineProblemToSolve>` has already been called with prefix **pr** in order to define the problem to solve.

In practice, this builds and performs the factorization of the local matrices used in the one level preconditioner.
The local matrices depend on the choice of :ref:`ffddmprecond <ffddmParametersGlobal>` and **VarfPrec**, see ``pr#aR``\ below.

**defines**:

-  ``string pr#prec`` equal to :ref:`ffddmprecond <ffddmParametersGlobal>`.
   Sets the type of one level preconditioner :math:`M^{-1}_1` to be used: “asm” (*Additive Schwarz*), “ras” (*Restricted Additive Schwarz*), “oras” (*Optimized Restricted Additive Schwarz*), “soras” (*Symmetric Optimized Restricted Additive Schwarz*) or “none” (no preconditioner).
-  ``matrix<pr#prfe#K>[int] pr#aR`` array (size ``prfe#prmesh#npart``) of local matrices used for the one level preconditioner.
   Each mpi rank of the spatial domain decomposition performs the :math:`LU` (or :math:`LDL^T`) factorization of the local matrix corresponding to its subdomain using the direct solver *MUMPS*.

   -  If **VarfPrec** is not a previously defined macro (just put *null* for example), the matrices ``pr#aR`` are set to be equal to the so-called local ‘Dirichlet’ matrices ``pr#aRd`` (see :ref:`ffddmsetupOperator <ffddmDocumentationDefineProblemToSolve>`).
      This is for the classical ASM preconditioner :math:`M^{-1}_1 = M^{-1}_{\text{ASM}} = \sum_{i=1}^N R_i^T A_i^{-1} R_i` or classical RAS preconditioner :math:`M^{-1}_1 = M^{-1}_{\text{RAS}} = \sum_{i=1}^N R_i^T D_i A_i^{-1} R_i` (it is assumed that :ref:`ffddmprecond <ffddmParametersGlobal>` is equal to “asm” or “ras”).
   -  If **VarfPrec** is a macro, it is assumed that **VarfPrec** defines an abstract bilinear form (see :ref:`ffddmsetupOperator <ffddmDocumentationDefineProblemToSolve>` for more details on how to define the abstract variational form as a macro).

      -  If :ref:`ffddmprecond <ffddmParametersGlobal>` is equal to “asm” or “ras”, the matrices ``pr#aR`` will be assembled as local ‘Dirichlet’ matrices in the same manner as ``pr#aRd``, but using the bilinear form defined by **VarfPrec** instead.
         This defines the ASM preconditioner as :math:`M^{-1}_1 = M^{-1}_{\text{ASM}} = \sum_{i=1}^N R_i^T {(A_i^{\text{Prec}})}^{-1} R_i` and the RAS preconditioner as :math:`M^{-1}_1 = M^{-1}_{\text{RAS}} = \sum_{i=1}^N R_i^T D_i {(A_i^{\text{Prec}})}^{-1} R_i`, where :math:`A_i^{\text{Prec}} = R_i A^{\text{Prec}} R_i^T`.
      -  If :ref:`ffddmprecond <ffddmParametersGlobal>` is equal to “oras” or “soras”, the matrices ``pr#aR`` will correspond to the discretization of the variational form **VarfPrec** in the subdomains :math:`\Omega_i`.
         In particular, various boundary conditions can be imposed at the interface between subdomains (corresponding to mesh boundary of label :ref:`ffddminterfacelabel <ffddmParametersGlobal>` set by the parent call to :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>`), such as Optimized Robin boundary conditions.
         We note the ORAS preconditioner as :math:`M^{-1}_1 = M^{-1}_{\text{ORAS}} = \sum_{i=1}^N R_i^T D_i {(B_i^{\text{Prec}})}^{-1} R_i` and the SORAS preconditioner as :math:`M^{-1}_1 = M^{-1}_{\text{SORAS}} = \sum_{i=1}^N R_i^T D_i {(B_i^{\text{Prec}})}^{-1} D_i R_i`.
-  ``func pr#prfe#K[int] pr#PREC1(pr#prfe#K[int] &ui)`` The function ``pr#PREC1`` computes the parallel application of the one level preconditioner :math:`M^{-1}_1`, i.e. the action of :math:`M^{-1}_1` on the local vector :math:`u_i`.
   In the sequential case, it computes the action of :math:`M^{-1}_1` on a global vector.
   The action of the inverse of local matrices ``pr#aRd`` is computed by forward-backward substitution using their :math:`LU` (or :math:`LDL^T`) decomposition.
-  ``func pr#prfe#K[int] pr#PREC(pr#prfe#K[int] &ui)`` The function ``pr#PREC`` corresponds to the action of the preconditioner :math:`M^{-1}` for problem **pr**.
   It coincides with the one level preconditioner ``pr#PREC1`` after the call to :ref:`ffddmsetupPrecond <ffddmDocumentationOneLevelPreconditioners>`.
   If a second level is subsequently added (see the next section about :ref:`Two level preconditioners <ffddmDocumentationTwoLevelPreconditioners>`), it will then coincide with the two level preconditioner :math:`M^{-1}_2` (see ``pr#PREC2level``).
-  ``func pr#prfe#K[int] pr#fGMRES(pr#prfe#K[int]& x0i, pr#prfe#K[int]& bi, real eps, int nbiter, string sprec)`` The function ``pr#fGMRES`` allows to solve the linear system :math:`A x = b` in parallel using the flexible GMRES method preconditioned by :math:`M^{-1}`.
   The action of the global operator :math:`A` is given by ``pr#A``, the action of the preconditioner :math:`M^{-1}` is given by ``pr#PREC`` and the scalar products are computed by ``pr#scalprod``.
   More details are given in the section :ref:`Solving the linear system <ffddmDocumentationSolvingLinearSystem>`.

.. raw:: html

   <!--
   ***For advanced users***:

   NewMacro pr#localmacroaug pr#prfe#prmesh#buildAug EndMacro
   IFMACRO(pr#localmacroaug,1)
   matrix<pr#prfe#K> pr#CSinterp;
   ENDIFMACRO
   -->

.. _ffddmDocumentationTwoLevelPreconditioners:

Two level preconditioners
-------------------------

The main ingredient of a two level preconditioner is the so-called ‘coarse space’ matrix :math:`Z`.

:math:`Z` is a rectangular matrix of size :math:`n \times n_c`, where usually :math:`n_c \ll n`.

:math:`Z` is used to build the ‘coarse space operator’ :math:`E = Z^T A Z`, a square matrix of size :math:`n_c \times n_c`.
We can then define the ‘coarse space correction operator’ :math:`Q = Z E^{-1} Z^T = Z (Z^T A Z)^{-1} Z^T`, which can then be used to enrich the one level preconditioner through a correction formula.
The simplest one is the *additive* coarse correction: :math:`M^{-1}_2 = M^{-1}_1 + Q`.
See ``pr#corr`` below for all other available correction formulas.

There are multiple ways to define a relevant coarse space :math:`Z` for different classes of problems.
:ref:`ffddmgeneosetup <ffddmDocumentationBuildingGeneoCoarseSpace>` defines a coarse space correction operator by building the GenEO coarse space, while :ref:`ffddmcoarsemeshsetup <ffddmDocumentationBuildingCoarseSpaceFromCoarseMesh>` builds the coarse space using a coarse mesh.

After a call to either :ref:`ffddmgeneosetup <ffddmDocumentationBuildingGeneoCoarseSpace>` or :ref:`ffddmcoarsemeshsetup <ffddmDocumentationBuildingCoarseSpaceFromCoarseMesh>`, the following variables and functions are set up:

-  ``int pr#ncoarsespace`` the size of the coarse space :math:`n_c`.
-  ``string pr#corr`` initialized with the value of :ref:`ffddmcorrection <ffddmParametersGlobal>`.
   Specifies the type of coarse correction formula to use for the two level preconditioner.
   The possible values are:

.. math::
    \begin{array}{llllll}
    \nonumber
        &&\text{"AD"}:&&\textit{Additive}, \quad &M^{-1} = M^{-1}_2 = \phantom{(I - Q A) }M^{-1}_1\phantom{ (I - A Q)} + Q\\
        &&\text{"BNN"}:&&\textit{Balancing Neumann-Neumann}, \quad &M^{-1} = M^{-1}_2 = (I - Q A) M^{-1}_1 (I - A Q) + Q\\
        &&\text{"ADEF1"}:&&\textit{Adapted Deflation Variant 1}, \quad &M^{-1} = M^{-1}_2 = \phantom{(I - Q A) }M^{-1}_1 (I - A Q) + Q\\
        &&\text{"ADEF2"}:&&\textit{Adapted Deflation Variant 2}, \quad &M^{-1} = M^{-1}_2 = (I - Q A) M^{-1}_1\phantom{ (I - A Q)} + Q\\
        &&\text{"RBNN1"}:&&\textit{Reduced Balancing Variant 1}, \quad &M^{-1} = M^{-1}_2 = (I - Q A) M^{-1}_1 (I - A Q)\\
        &&\text{"RBNN2"}:&&\textit{Reduced Balancing Variant 2}, \quad &M^{-1} = M^{-1}_2 = (I - Q A) M^{-1}_1\phantom{ (I - A Q)}\\
        &&\text{"none"}:&&\textit{no coarse correction}, \quad &M^{-1} = M^{-1}_2 = \phantom{(I - Q A) }M^{-1}_1\phantom{ (I - A Q)}\\
    \end{array}

-  Note that *AD*, *ADEF1* and *RBNN2* only require one application of :math:`Q`, while *BNN*, *ADEF2* and *RBNN1* require two.
   The default coarse correction is *ADEF1*, which is cheaper and generally as robust as *BNN* or *ADEF2*.
-  ``func pr#prfe#K[int] pr#Q(pr#prfe#K[int] &ui)`` The function ``pr#Q`` computes the parallel application of the coarse correction operator :math:`Q`, i.e. the action of :math:`Q = Z E^{-1} Z^T` on the local vector :math:`u_i`.
   In the sequential case, it computes the action of :math:`Q` on a global vector.
   The implementation differs depending on the method used to build the coarse space (with GenEO or using a coarse mesh), but the idea is the same: the action of the transpose of the distributed operator :math:`Z` on the distributed vector :math:`u_i` is computed in parallel, with the contribution of all subdomains being gathered in a vector of size :math:`n_c` in the mpi process of rank 0.
   The action of the inverse of the coarse space operator :math:`E` is then computed by forward-backward substitution using its :math:`LU` (or :math:`LDL^T`) decomposition previously computed by the first ``pr#prfe#prmesh#pCS`` ranks of the mpi communicator.
   The result is then sent back to all subdomains to perform the last application of :math:`Z` and obtain the resulting local vector in each subdomain.
-  ``func pr#prfe#K[int] pr#PREC2level(pr#prfe#K[int] &ui)`` The function ``pr#PREC2level`` computes the parallel application of the two level preconditioner :math:`M^{-1}_2`, i.e. the action of :math:`M^{-1}_2` on the local vector :math:`u_i`.
   In the sequential case, it computes the action of :math:`M^{-1}_2` on a global vector.
   The two level preconditioner depends on the choice of the coarse correction formula which is determined by ``pr#corr``, see above.

.. raw:: html

   <!--
   ***For advanced users***:

   int pr#bCM = 0;
   -->

.. _ffddmDocumentationBuildingGeneoCoarseSpace:

Building the GenEO coarse space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: freefem
   :linenos:

   ffddmgeneosetup(pr,Varf)

This builds the GenEO coarse space for problem **pr**.
This will create and expose variables whose names will be prefixed by **pr**, see below.
It is assumed that :ref:`ffddmsetupPrecond <ffddmDocumentationOneLevelPreconditioners>` has already been called for prefix **pr** in order to define the one level preconditioner for problem **pr**.
The GenEO coarse space is :math:`Z = (R_i^T D_i V_{i,k})^{i=1,...,N}_{\lambda_{i,k} \ge \tau}`, where :math:`V_{i,k}` are eigenvectors corresponding to eigenvalues :math:`\lambda_{i,k}` of the following local generalized eigenvalue problem in subdomain :math:`i`:

:math:`D_i A_i D_i V_{i,k} = \lambda_{i,k} A_i^{\text{Neu}} V_{i,k}`,

where :math:`A_i^{\text{Neu}}` is the local Neumann matrix of subdomain :math:`i` (with Neumann boundary conditions at the subdomain interface).

In practice, this builds and factorizes the local Neumann matrices :math:`A_i^{\text{Neu}}` corresponding to the abstract bilinear form given by the macro **Varf** (see :ref:`ffddmsetupOperator <ffddmDocumentationDefineProblemToSolve>` for more details on how to define the abstract variational form as a macro).
In the GenEO method, the abstract bilinear form **Varf** is assumed to be the same as the one used to define the problem **pr** through the previous call to :ref:`ffddmsetupOperator <ffddmDocumentationDefineProblemToSolve>`.
The local generalized eigenvalue problem is then solved in each subdomain to find the eigenvectors :math:`V_{i,k}` corresponding to the largest eigenvalues :math:`\lambda_{i,k}` (see ``pr#Z`` below).
The number of computed eigenvectors :math:`\nu` is given by :ref:`ffddmnu <ffddmParametersGlobal>`.
The eigenvectors selected to enter :math:`Z` correspond to eigenvalues :math:`\lambda_{i,k}` larger than :math:`\tau`, where the threshold parameter :math:`\tau` is given by :ref:`ffddmtau <ffddmParametersGlobal>`.
If :ref:`ffddmtau <ffddmParametersGlobal>` :math:`= 0`, all :ref:`ffddmnu <ffddmParametersGlobal>` eigenvectors are selected.
Finally, the coarse space operator :math:`E = Z^T A Z` is assembled and factorized (see ``pr#E`` below).

**defines**:

-  ``pr#prfe#K[int][int] pr#Z`` array of local eigenvectors :math:`Z_{i,k} = D_i V_{i,k}` obtained by solving the local generalized eigenvalue problem above in the subdomain of this mpi rank using *Arpack*.
   The number of computed eigenvectors :math:`\nu` is given by :ref:`ffddmnu <ffddmParametersGlobal>`.
   The eigenvectors selected to enter :math:`Z` correspond to eigenvalues :math:`\lambda_{i,k}` larger than :math:`\tau`, where the threshold parameter :math:`\tau` is given by :ref:`ffddmtau <ffddmParametersGlobal>`.
   If :ref:`ffddmtau <ffddmParametersGlobal>` :math:`= 0`, all :ref:`ffddmnu <ffddmParametersGlobal>` eigenvectors are selected.
-  ``matrix<pr#prfe#K> pr#E`` the coarse space operator :math:`E = Z^T A Z`.
   The matrix ``pr#E`` is assembled in parallel and is factorized by the parallel direct solver *MUMPS* using the first ``pr#prfe#prmesh#pCS`` ranks of the mpi communicator, with mpi rank 0 as the master process.
   The number of mpi processes dedicated to the coarse problem is set by the underlying mesh decomposition of problem **pr**, which also specifies if these mpi ranks are excluded from the spatial decomposition or not.
   These parameters are set by :ref:`ffddmpCS <ffddmParametersGlobal>` and :ref:`ffddmexclude <ffddmParametersGlobal>` when calling :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>` (see :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>` for more details).

.. raw:: html

   <!--
   ***For advanced users***:

   int pr#si;

   pr#sizelg(pr#prfe#prmesh#npart), pr#offseti(pr#prfe#prmesh#npart);

   int[int] pr#sizelgworld(mpiSize(pr#prfe#prmesh#mpicomm)), pr#offsetiworld(mpiSize(pr#prfe#prmesh#mpicomm));

   matrix<pr#prfe#K> pr#matN;
   -->

.. _ffddmDocumentationBuildingCoarseSpaceFromCoarseMesh:

Building the coarse space from a coarse mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: freefem
   :linenos:

   ffddmcoarsemeshsetup(pr,Thc,VarfEprec,VarfAprec)

builds the coarse space for problem **pr** from a coarse problem which corresponds to the discretization of a variational form on a coarser mesh **Thc** of :math:`\Omega`.
This will create and expose variables whose names will be prefixed by **pr**, see below.
It is assumed that :ref:`ffddmsetupPrecond <ffddmDocumentationOneLevelPreconditioners>` has already been called for prefix **pr** in order to define the one level preconditioner for problem **pr**.
The abstract variational form for the coarse problem can differ from the original problem **pr** and is given by macro **VarfEprec** (see :ref:`ffddmsetupOperator <ffddmDocumentationDefineProblemToSolve>` for more details on how to define the abstract variational form as a macro).
For example, absorption can be added in the preconditioner for wave propagation problems, see examples for Helmholtz and Maxwell equations in the :ref:`Examples <ffddmExamples>` section.

The coarse space :math:`Z` corresponds to the interpolation operator from the coarse finite element space to the original finite element space of the problem.
Thus, the coarse space operator :math:`E = Z^T A^{\text{Eprec}} Z` corresponds to the matrix of the problem given by **VarfEprec** discretized on the coarse mesh **Thc** and is assembled as such.

Similarly, **VarfAprec** specifies the global operator involved in multiplicative coarse correction formulas.
For example, :math:`M^{-1}_{2,\text{ADEF1}} = M^{-1}_1 (I - A^{\text{Aprec}} Q) + Q` (where :math:`Q = Z E^{-1} Z^T`).
:math:`A^{\text{Aprec}}` defaults to :math:`A` if **VarfAprec** is not a valid macro (you can put *null* for example).

**defines**:

-  ``meshN pr#ThCoarse`` the coarse mesh **Thc**
-  ``fespace pr#VhCoarse`` the coarse finite element space of type ``pr#prfe#fPk`` defined on the coarse mesh ``pr#ThCoarse``
-  ``matrix<pr#prfe#K> pr#AglobEprec`` the global matrix :math:`A^{\text{Aprec}}` corresponding to the discretization of the variational form given by the macro **VarfAprec** on the global finite element space ``pr#prfe#Vhglob``.
   Defined only in the sequential case.
   ``pr#AglobEprec`` is equal to ``pr#Aglobal`` if **VarfAprec** is not a valid macro.
-  ``matrix<pr#prfe#K> pr#aRdEprec`` the local ‘Dirichlet’ matrix corresponding to **VarfAprec**; it is the local restriction of the global operator :math:`A^{\text{Aprec}}` to the subdomain, equivalent to :math:`A^{\text{Aprec}}_i = R_i A^{\text{Aprec}} R_i^T` with :math:`A^{\text{Aprec}}` the global matrix corresponding to the discretization of the variational form given by the macro **VarfAprec** on the global finite element space.
   Defined only if this mpi rank is not excluded from the spatial domain decomposition, i. e. ``prmesh#excluded`` = 0.
   ``pr#aRdEprec`` is equal to ``pr#aRd[mpiRank(prmesh#commddm)]`` if **VarfAprec** is not a valid macro.
-  ``func pr#prfe#K[int] pr#AEprec(pr#prfe#K[int] &ui)`` The function ``pr#AEprec`` computes the parallel matrix-vector product, i.e. the action of the global operator :math:`A^{\text{Aprec}}` on the local vector :math:`u_i`.
   The computation is equivalent to :math:`R_i (\sum_{j=1}^N R_j^T D_j A^{\text{Aprec}}_j u_j)` and is performed in parallel using local matrices ``pr#aRdEprec`` and the function ``pr#prfe#update``.
   In the sequential case, the global matrix ``pr#AglobEprec`` is used instead.
-  ``matrix<pr#prfe#K> pr#ZCM`` the interpolation operator :math:`Z` from the coarse finite element space ``pr#VhCoarse`` to the global finite element space ``pr#prfe#Vhglob``.
   Defined only in the sequential case.
-  ``matrix<pr#prfe#K> pr#ZCMi`` the local interpolation operator :math:`Z_i` from the coarse finite element space ``pr#VhCoarse`` to the local finite element space ``pr#prfe#Vhi``.
   Defined only if this mpi rank is not excluded from the spatial domain decomposition, i. e. ``prmesh#excluded`` = 0.
   ``pr#ZCMi`` is used for the parallel application of :math:`Z` and :math:`Z^T`.
-  ``matrix<pr#prfe#K> pr#ECM`` the coarse space operator :math:`E = Z^T A^{\text{Eprec}} Z`.
   The matrix ``pr#ECM`` is assembled by discretizing the variational form given by **VarfEprec** on the coarse mesh and factorized by the parallel direct solver *MUMPS* using the first ``pr#prfe#prmesh#pCS`` ranks of the mpi communicator, with mpi rank 0 as the master process.
   The number of mpi processes dedicated to the coarse problem is set by the underlying mesh decomposition of problem **pr**, which also specifies if these mpi ranks are excluded from the spatial decomposition or not.
   These parameters are set by :ref:`ffddmpCS <ffddmParametersGlobal>` and :ref:`ffddmexclude <ffddmParametersGlobal>` when calling :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>` (see :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>` for more details).

.. _ffddmDocumentationSolvingLinearSystem:

Solving the linear system
-------------------------

.. code-block:: freefem
   :linenos:

   func pr#prfe#K[int] pr#fGMRES(pr#prfe#K[int]& x0i, pr#prfe#K[int]& bi, real eps, int itmax, string sp)

solves the linear system for problem **pr** using the flexible GMRES algorithm with preconditioner :math:`M^{-1}` (corresponding to ``pr#PREC``).
Returns the local vector corresponding to the restriction of the solution to ``pr#prfe#Vhi``.
**x0i** and **bi** are local distributed vectors corresponding respectively to the initial guess and the right-hand side (see :ref:`ffddmbuildrhs <ffddmDocumentationBuildRhs>`).
**eps** is the stopping criterion in terms of the relative decrease in residual norm.
If **eps** :math:`< 0`, the residual norm itself is used instead.
**itmax** sets the maximum number of iterations.
**sp** selects between the ``"left"`` or ``"right"`` preconditioning variants: *left* preconditioned GMRES solves :math:`M^{-1} A x = M^{-1} b`, while *right* preconditioned GMRES solves :math:`A M^{-1} y = b` for :math:`y`, with :math:`x = M^{-1} y`.

.. _ffddmDocumentationHPDDMffddm:

Using *HPDDM* within *ffddm*
----------------------------

**ffddm** allows you to use **HPDDM** to solve your problem, effectively replacing the **ffddm** implementation of all parallel linear algebra computations.
**ffddm** can then be viewed as a finite element interface for **HPDDM**.

You can use **HPDDM** features unavailable in **ffddm** such as advanced Krylov subspace methods implementing block and recycling techniques.

To switch to **HPDDM**, simply define the macro ``pr#withhpddm`` before using :ref:`ffddmsetupOperator <ffddmDocumentationDefineProblemToSolve>`. You can then pass **HPDDM** options
with command-line arguments or directly to the underlying **HPDDM** operator ``pr#hpddmOP``. Options need to be prefixed by the operator prefix:

.. code-block:: freefem
  :linenos:

  macro PBwithhpddm()1 // EOM
  ffddmsetupOperator( PB , FE , Varf )
  set(PBhpddmOP,sparams="-hpddm_PB_krylov_method gcrodr");

You can also choose to replace only the Krylov solver, by defining the macro ``pr#withhpddmkrylov`` before using :ref:`ffddmsetupOperator <ffddmDocumentationDefineProblemToSolve>`.
Doing so, a call to ``pr#fGMRES`` will call the **HPDDM** Krylov solver, with **ffddm** providing the operator and preconditioner through ``pr#A`` and ``pr#PREC``.

An example can be found in **Helmholtz-2d-HPDDM-BGMRES.edp**, see the :ref:`Examples <ffddmExamples>` section.

.. _ffddmDocumentationAdvanced:

Advanced use
----------------------------

.. raw:: html

  <!--
  .. _ffddmDocumentationNonlinearTimedependent:

  Nonlinear and time dependent problems
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  .. code-block:: freefem
    :linenos:

    pr#fromVhi(ui,VhName,res)

  .. code-block:: freefem
    :linenos:

  ffddmbuildDmeshAug(pr,Th,comm)
  -->


.. _ffddmDocumentationPartitionUnityEdge:

Local finite element spaces for non Lagrange finite elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For Lagrange finite elements, the partition of unity :math:`(D_i)_{i=1,...,N}` (see ``pr#Dk`` and ``pr#Dih``) is built by interpolating the local P1 partition of unity function onto the components of the **Pk** finite element space ``pr#Vhi``.
For non Lagrange finite element spaces, such as Raviart–Thomas or Nédélec edge elements, the definition of the degrees of freedom can be more involved, and interpolating the P1 partition of unity functions directly is inappropriate.
The idea is then to use a "pseudo" finite element **Pkpart** derived from **Pk** which is suitable for interpolating the P1 partition of unity, in the sense that it will produce a partition of unity for **Pk**.

For example, for first-order Nédélec edge elements (*Edge03d*), whose degrees of freedom are the circulations along the edges, we define the "pseudo" finite element *Edge03ds0* which can be seen as a scalar Lagrange counterpart: the numbering of the degrees of freedom is the same, but they correspond to the value at the edge midpoints.

For Lagrange finite elements, the distributed finite element spaces are built using :ref:`ffddmbuildDfespace <ffddmDocumentationLocalFiniteElementSpaces>`. Here you must use **ffddmbuildDfespaceEdge**, which builds the distributed finite element space using a "pseudo" finite element to build the partition of unity:

.. code-block:: freefem
   :linenos:

   ffddmbuildDfespaceEdge(pr,prmesh,scalar,def,init,Pk,defpart,initpart,Pkpart)

where macros **defpart** and **initpart** specify how to define and interpolate a function in the 'pseudo' finite element space **Pkpart**, similar to **def** and **init** for **Pk**.

An example with first-order Nédélec edge elements (*Edge03d* + *Edge03ds0*) for Maxwell equations can be found in **Maxwell-3d-simple.edp**, see the :ref:`Examples <ffddmExamples>` section.

.. _ffddmDocumentationInexactCoarseSolve:

Inexact coarse solves for two level methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have seen in the :ref:`Two level preconditioners section <ffddmDocumentationTwoLevelPreconditioners>` that two level methods produce a ‘coarse space operator’ :math:`E` that needs to be inverted at each iteration.
By default the coarse space operator matrix is factorized by the direct solver *MUMPS*. This can become a bottleneck and hinder scalability for large problems, where :math:`E` can become too large to be factorized efficiently.
To remedy this, we can instead opt to use an iterative method to solve the coarse problem at each iteration.  Moreover, in order to retain robustness, a DD preconditioner can be used to solve the inner coarse problem more efficiently.

.. raw:: html

  <!--
  Three level GenEO
  '''''''''''''''''
  -->

Coarse mesh and inexact coarse solve
''''''''''''''''''''''''''''''''''''''''

When the coarse problem comes from a coarse mesh discretization, a natural way to do inexact coarse solve is to use a one level domain decomposition method on the coarse problem, with the same subdomain partitioning for the coarse and fine meshes.
This means that each processor is associated to one spatial subdomain and hosts the two local (nested) coarse and fine submeshes corresponding to this subdomain, as well as the corresponding local matrices for the two discretizations.
This natural choice offers interesting benefits: 

-  We naturally recover a load-balanced parallel implementation, provided that the initial partitioning is balanced.
-  The communication pattern between neighboring subdomains is the same for the coarse and fine discretizations.
-  The assembly and the application of the interpolation operator :math:`Z` (and :math:`Z^T`) between the fine and the coarse spaces can be computed locally in each subdomain and require no communication.

In **ffddm**, the first step is to build the two nested mesh decompositions using **ffddmbuildDmeshNested**:

.. code-block:: freefem
   :linenos:

   ffddmbuildDmeshNested(pr,Thc,s,comm)

decomposes the coarse mesh **Thc** into overlapping submeshes and creates the fine decomposition by locally refining submeshes by a factor of **s**, i.e. splitting each mesh element into :math:`s^d` elements, :math:`s \geq 1`.
This will create and expose variables corresponding to both decompositions, prefixed by **pr** for the fine mesh and by **pr#Coarse** for the coarse mesh (see :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>`).  
It also sets the integer variable ``pr#binexactCS`` to 1, which specifies that any two level method defined on mesh prefix **pr** will use inexact coarse solves.

The distributed finite element spaces, operators and preconditioners can then be defined for both decompositions. Here is an example where the coarse problem is solved using a one level method:

.. code-block:: freefem
   :linenos:

   ffddmbuildDmeshNested(M, Thc, 3, mpiCommWorld)

   ffddmbuildDfespace(FE, M, real, def, init, Pk)
   ffddmbuildDfespace(FECoarse, MCoarse, real, def, init, Pk)

   // coarse operator (Varf of E):
   ffddmsetupOperator(PBCoarse, FECoarse, VarfEprec)
   // one level preconditioner for the coarse problem:
   ffddmsetupPrecond(PBCoarse, VarfPrecC)

   // operator for the fine problem:
   ffddmsetupOperator(PB, FE, Varf)
   // one level preconditioner for the fine problem:
   ffddmsetupPrecond(PB, VarfPrec)

   // add the second level:
   ffddmcoarsemeshsetup(PB, Thc, VarfEprec, null)

   [...]
   u[] = PBfGMRES(x0, rhs, 1.e-6, 200, "right");

**Remarks**:

- Note that the different prefixes need to match: prefixes for the coarse decomposition have to be those of the fine decomposition, appended with ``Coarse``.
- The operator and preconditioner for the coarse problem have to be defined before those of the fine problem, because the ``pr#Q`` function is actually defined by ``ffddmsetupPrecond`` and involves a call to ``pr#CoarsefGMRES`` (which is defined by ``ffddmsetupPrecond`` for the coarse problem) for the iterative solution of the coarse problem if ``pr#prfe#prmesh#binexactCS`` :math:`\neq 0`.
- In this case, ``ffddmcoarsemeshsetup`` does not use **Thc** or **VarfEprec** and only builds the local interpolation matrices between fine and coarse local finite element spaces ``pr#prfe#Vhi`` and ``pr#prfe#CoarseVhi`` to be able to apply :math:`Z` and :math:`Z^T`.
- The GMRES tolerance for the inner solution of the coarse problem is set by :ref:`ffddminexactCStol <ffddmParametersGlobal>` and is equal to 0.1 by default.

In practice, these methods can give good results for wave propagation problems, where the addition of artificial absorption in the preconditioner helps with the convergence of the one level method for the inner solution of the coarse problem.
You can find an example for Maxwell equations in **Maxwell_Cobracavity.edp**, see the :ref:`Examples <ffddmExamples>` section. More details can be found `here`_ and in 

  \M. Bonazzoli, V. Dolean, I. G. Graham, E. A. Spence, P.-H. Tournier. Domain decomposition preconditioning for the high-frequency time-harmonic Maxwell equations with absorption. Mathematics of Computation, 2019. DOI: https://doi.org/10.1090/mcom/3447 

.. _here: ../../_static/html/tutorial-slides.html#26

.. raw:: html

  <!--
  Computing integrals
  ~~~~~~~~~~~~~~~~~~~
  -->



