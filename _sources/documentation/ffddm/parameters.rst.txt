Parameters
==========

.. _ffddmParametersCommandLine:

Command-line arguments
----------------------

-  :freefem:`-ffddm_verbosity N`, the level of verbosity of **ffddm**, see :ref:`ffddmverbosity <ffddmParametersGlobal>` (default 3).
-  :freefem:`-seqddm N` use **ffddm** in sequential mode, with N the number of subdomains.
-  :freefem:`-noGlob` if present, do not define any global quantity (such as saving the global mesh for plotting or building the global restriction matrices).
   Cannot be used in sequential mode or with plotting.
-  :freefem:`-ffddm_partitioner N` specifies how to partition the initial domain, see :ref:`ffddmpartitioner <ffddmParametersGlobal>` (default 1, *metis*).
-  :freefem:`-ffddm_overlap N` specifies the width of the overlap region between subdomains, see :ref:`ffddmoverlap <ffddmParametersGlobal>` (default 1).
-  :freefem:`-ffddm_master_p N`, number of master processes for the coarse problem (for two level preconditioners), see :ref:`ffddmpCS <ffddmParametersGlobal>` (default 1).
-  :freefem:`-ffddm_master_exclude 0|1` exclude master processes from the domain decomposition, see :ref:`ffddmexclude <ffddmParametersGlobal>` (default 0).
-  :freefem:`-ffddm_split N`, level of refinement of the local submeshes with respect to the initial global mesh, see :ref:`ffddmsplit <ffddmParametersGlobal>` (default 1).
-  :freefem:`-ffddm_schwarz_method S`, specifies the type of one level preconditioner :math:`M^{-1}_1`: “asm” (*Additive Schwarz*), “ras” (*Restricted Additive Schwarz*), “oras” (*Optimized Restricted Additive Schwarz*), “soras” (*Symmetric Optimized Restricted Additive Schwarz*) or “none” (no preconditioner), see :ref:`ffddmprecond <ffddmParametersGlobal>` (default “ras”).
-  :freefem:`-ffddm_geneo_nu N`, number of local eigenvectors to compute in each subdomain when solving the local generalized eigenvalue problem for the GenEO method, see :ref:`ffddmnu <ffddmParametersGlobal>` (default 20).
-  :freefem:`-ffddm_geneo_threshold R`, threshold parameter for selecting local eigenvectors when solving the local generalized eigenvalue problems for the GenEO method, see :ref:`ffddmtau <ffddmParametersGlobal>` (default 0.5).
   If the command-line parameter **-ffddm_geneo_nu N** is used, then :ref:`ffddmtau <ffddmParametersGlobal>` is initialized to 0.
-  :freefem:`-ffddm_schwarz_coarse_correction S`, specifies the coarse correction formula to use for the two level preconditioner: “AD” (*Additive*), “BNN” (*Balancing Neumann-Neumann*), “ADEF1” (*Adapted Deflation Variant 1*), “ADEF2” (*Adapted Deflation Variant 2*), “RBNN1” (*Reduced Balancing Variant 1*), “RBNN2” (*Reduced Balancing Variant 2*) or “none” (no coarse correction), see :ref:`ffddmcorrection <ffddmParametersGlobal>` (default “ADEF1”).
-  :freefem:`-ffddm_inexactCS_tol R`, specifies the GMRES tolerance for the inner solution of the coarse problem when using a two level method with approximate coarse solves, see :ref:`ffddminexactCStol <ffddmParametersGlobal>` (default 0.1).

.. _ffddmParametersGlobal:

Global parameters
-----------------

-  :freefem:`ffddmverbosity` initialized by command-line argument **-ffddm_verbosity N**, specifies the level of verbosity of **ffddm** (default 3).
-  :freefem:`ffddmpartitioner` initialized by command-line argument **-ffddm_partitioner N**, specifies how to partition the initial domain:

   -  N=0: user-defined partition through the definition of a macro, see :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>`
   -  N=1: use the automatic graph partitioner *metis* (default)
   -  N=2: use the automatic graph partitioner *scotch*
-  :freefem:`ffddmoverlap` initialized by command-line argument **-ffddm_overlap N**, specifies the number of layers of mesh elements in the overlap region between subdomains ; N >= 0 (default 1).
   **Remark** The actual width of the overlap region between subdomains is 2N, since each subdomain is extended by N layers of elements in a symmetric way.  **Remark 2** if :ref:`ffddmoverlap <ffddmParametersGlobal>` = 0, the construction is a bit different, since only interface unknowns are shared. In that case, the user is required to define the macro :freefem:`prmesh#mminoverlap` to 1 before calling :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>`. It can be done e.g. on the command line: :freefem:`-ffddm_overlap 0 -DMESHmminoverlap=1`, where :freefem:`MESH` is your distributed mesh prefix.
-  :freefem:`ffddminterfacelabel` the label of the new border of the subdomain meshes (the interface between the subdomains) (default 10).
   Used for imposing problem-dependent boundary conditions at the interface between subdomains for the preconditioner, for example optimized Robin boundary conditions (see ORAS).
-  :freefem:`ffddmpCS` initialized by command-line argument **-ffddm_master_p N**, number of mpi processes used for the assembly and resolution of the coarse problem for two level preconditioners (default 1).
-  :freefem:`ffddmexclude` initialized by command-line argument **-ffddm_master_exclude**, 0 or 1 (default 0).
   If true, mpi ranks participating in the assembly and resolution of the coarse problem for two level preconditioners will be excluded from the spatial domain decomposition and will only work on the coarse problem.
-  :freefem:`ffddmsplit` initialized by command-line argument **ffddm_split N**, level of refinement of the local submeshes with respect to the initial global mesh (default 1).
   This is useful for large problems, where we want to avoid working with a very large global mesh.
   The idea is to start from a coarser global mesh, and generate finer local meshes in parallel during the mesh decomposition step in order to reach the desired level of refinement for the subdomains.
   For example, calling :ref:`ffddmbuildDmesh <ffddmDocumentationOverlappingMeshDecomposition>` with :ref:`ffddmsplit <ffddmParametersGlobal>` = 3 will generate local submeshes where each mesh element of the initial mesh is split into :math:`3^d` elements.
-  :freefem:`ffddmprecond` initialized by command-line argument **-ffddm_schwarz_method S**, specifies the type of one level preconditioner :math:`M^{-1}_1` to build when calling :ref:`ffddmsetupPrecond <ffddmDocumentationOneLevelPreconditioners>`: “asm” (*Additive Schwarz*), “ras” (*Restricted Additive Schwarz*), “oras” (*Optimized Restricted Additive Schwarz*), “soras” (*Symmetric Optimized Restricted Additive Schwarz*) or “none” (no preconditioner).
   Default is “ras”.
   See :ref:`ffddmsetupPrecond <ffddmDocumentationOneLevelPreconditioners>` for more details.
-  :freefem:`ffddmnu` initialized by command-line argument **-ffddm_geneo_nu N**, number of local eigenvectors to compute in each subdomain when solving the local generalized eigenvalue problem for the GenEO method (default 20).
   See :ref:`ffddmgeneosetup <ffddmDocumentationBuildingGeneoCoarseSpace>` for more details.
-  :freefem:`ffddmtau` initialized by command-line argument **-ffddm_geneo_threshold R**, threshold parameter for selecting local eigenvectors when solving the local generalized eigenvalue problems for the GenEO method (default 0.5).
   If the command-line parameter **-ffddm_geneo_nu N** is used, then :ref:`ffddmtau <ffddmParametersGlobal>` is initialized to 0.
   See :ref:`ffddmgeneosetup <ffddmDocumentationBuildingGeneoCoarseSpace>` for more details.
-  :freefem:`ffddmcorrection` initialized by command-line argument **-ffddm_schwarz_coarse_correction S**, specifies the coarse correction formula to use for the two level preconditioner: “AD” (*Additive*), “BNN” (*Balancing Neumann-Neumann*), “ADEF1” (*Adapted Deflation Variant 1*), “ADEF2” (*Adapted Deflation Variant 2*), “RBNN1” (*Reduced Balancing Variant 1*), “RBNN2” (*Reduced Balancing Variant 2*) or “none” (no coarse correction).
   Default is “ADEF1”.
   See the section about :ref:`Two level preconditioners <ffddmDocumentationTwoLevelPreconditioners>` for more details.
-  :freefem:`ffddminexactCStol` initialized by command-line argument **-ffddm_inexactCS_tol R**, GMRES tolerance for the inner solution of the coarse problem when using a two level method with approximate coarse solves (default 0.1).
   See the section about :ref:`Approximate coarse solves for two level methods <ffddmDocumentationApproximateCoarseSolve>` for more details.
