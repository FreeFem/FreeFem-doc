.. role:: freefem(code)
   :language: freefem

.. _new-features:

New features
============

The notable changes of each **FreeFEM** release are listed below.

Version 4.15 (6 December 2024)
------------------------------

* Added

  - FreeFEM can now run Markdown (.md) files as well as .edp files. Markdown can be used to document your FreeFEM scripts (see for example the .md files in the `examples/examples <https://github.com/FreeFem/FreeFem-sources/tree/master/examples/examples>`__ directory). When running a .md file, FreeFEM will execute the code contained in the FreeFEM code blocks, delimited by

    .. code-block:: markdown
      :linenos:

      ```freefem
      [freefem code]
      ```

    Documented Markdown examples in the distribution can be viewed on the `documentation website <https://doc.freefem.org>`__ (toggle Search in examples). Markdown can be previewed with e.g Visual Studio Code, with the `FreeFEM VS Code extension <https://marketplace.visualstudio.com/items?itemName=Pierre-Marchand.vscode-freefem>`__ providing syntax highlighting for FreeFEM code blocks.

  - add command :console:`md2edp` to extract freefem code form Markdown file.

  - add read binary file with short number in bfstream plugin
  - add 3d case in ClosePoints plugin : function Voisinage3, ClosePoints3

  - add :freefem:`real[int,int] at=A';`

  - Powell Sabin :freefem:`splitmesh6PowellSabin(Th)`  splitting for Scott–Vogelius lowest Stokes Element in 2D
    in plugin/seq/splitmesh6.cpp

  - Worsey Farin :freefem:`splitmesh12WorseyFarin(Th3)` splitting for Scott–Vogelius lowest Stokes Element in 3D
    in  plugin/seq/splitmesh12.cpp

  - functional interface with :freefem:`fgmres` (Linear and affine) in real and complex case
    see tutorial/algo.edp

* Changed

  - PETSc 3.22.2
  - move the plugin msh3 in kernel , so remove all :freefem:`load "msh3"` in all examples and .idp files

* Fixed

  - try to fix orientation of internal edges in 2d
  - fix missing break and continue in :freefem:`for [i,ai:a]` loop
  - line number in macro error (thanks to P-H Tournier)
  - fix problem in integration of moving test or unknown function in 2d mesh:

    .. code-block:: freefem
      :linenos:

      varf ab([u],[v]) = int2d(Th,mapu=[Xo,Yo])(u*v);
      matrix AB = ab(Zh,Rh);

    we remove a piece of code.
  - fix problem in mesh of a ring from a square (missing option to remove duplicate vertices)
    example: :freefem:`mesh Th2=square(19,5,[(1+y)*cos(2*pi*x),(1+y)*sin(x*pi*2)],removeduplicate=1);`
  - Segmentation fault problem in Write_hdf5 (problem of allocation in the stack not on the heap) if the mesh2 is too large
  - fix :freefem:`(A'+A)` where func  :freefem:`A =[[0,1],[0,0]];`
  - correct integer overflow (in rare case) when calling INTER_SEG1d use in interpolation on Th3,ThS,ThL
    thanks to G. Sadaka for finding the bug
  - correct hidden faces on surface mesh (ffglut)
  - remove optimisation flag ppm2rnm.cpp in macOS (load trap)
  - correct abcisse curviline on reparametrage function (setcurveabcisse) in :freefem:`load "Curvature"`
  - correct compilation on ARM sonoma xcode 15.2
  - correct mpi essai.edp MPI_ANY_SOURCE Pb

Version 4.14
------------

* Added

  - Finite element BDM2 and BDM2ortho in test, Bug in BDM2ortho corrected the 4 sept 2014 in version: v4.13-130-g1af52457
  - Conversion of matrix or transpose of matrix in :freefem:`int[int][int]` array to get the structure of sparse matrix.
    see tutorial/sparse-matrix.edp example at end

    .. code-block:: freefem
      :linenos:

      matrix A = va(Ph,Vh);
      int[int] a = A, at= A';

  - a meshL finite function can be see as real function with 1, or 2 parameters

    .. code-block:: freefem
      :linenos:

      meshL ThL = segment(10); fespace VhL(ThL,P1); VhL u= x;
      cout << u(0.5)   << endl;
      cout << u(0.5,0) << endl;

  - Exemple to code convolution of 2 functions with one with a small support to be not too expensive
    see tutorial/Convolution-Sample.edp example
  - Support for dense blocks in PETSc matrices
  - GenEO for saddle-point examples with PCHPDDM in PETSc
  - Distributed ParaView output on :freefem:`meshS`
  - Interface to :freefem:`mmg2d` for two-dimensional :freefem:`mesh`
  - Support for Mmg parameters :freefem:`localParameter`, :freefem:`-nosizreq`, :freefem:`-hgradreq`

* Changed

  - PETSc 3.20.2

* Fixed

  - bug in P3pnc3d in vectorial case (thanks to loic.balaziatchynillama@cea.fr)
  - in segment(10,region=1,label=ll); region is now used..

Version 4.13
------------

* Added

  - Composite FE spaces and variational forms for coupled problems (see :ref:`Composite finite element spaces <composite>` ):

    - can now define composite FE spaces with different meshes/mesh types as

      .. code-block:: freefem
        :linenos:

        fespace Uh(Th1,[P2,P2]);
        fespace Ph(Th2,P1);
        fespace Vh=Uh*Ph;

    - can define coupled problems using composite FE spaces, or directly with `<` `>` syntax:

      .. code-block:: freefem
        :linenos:

        fespace Uh(Th1,[P2,P2]);
        fespace Ph(Th2,P1);
        Uh [u1,u2],[v1,v2];
        Ph p,q;

        solve Pb (<[u1,u2],[p]>, <[v1,v2],[q]>) = ...

      see `examples/examples/stokes_composite.edp` and `examples/examples/stokes_periodic_composite.edp`
    - this new type of composite problem can be used for FEM-BEM coupling and also benefits from automatic parallel assembly (in test) and can be easily solved using the distributed solver MUMPS, see `examples/bem/Helmholtz-2d-FEM-BEM-coupling-MUMPS-composite.edp`
    - composite problems can also be solved using PETSc (in test), see `examples/hpddm/Helmholtz-2d-FEM-BEM-coupling-PETSc-composite.edp`

  - remove spurious cout in Curve/Line DG definition.

  - | add New Finite element 2d on mesh :  RT0dc (discontinuous RT0 ) in plugin Element_Mixte
    | see example plugin/RT0dc.edp
    | and P1nc (Crouziex-Raviat) + bulle : name P1bnc in plugin Element_P1ncdc
    | and P1nc totally discontinous + bulle  ; name P1bdcnc in plugin Element_P1ncdc
    | see example plugin/example testp1dcnc.edp
    | for akram.beni-hamad@inria.fr

  - | add New finite element:  P4S P4 on meshS , P3pnc3d in Element_P3pnc_3d (Couziex-Raviart with P3 )
    | see loic.balaziatchynillama@cea.fr for more information.

  - | add new interface for metis (see examples/plugin/metis.edp)
  - | Correct jump, mean, otherside of finite element function on mesh3, meshS, meshL
    | (add missing code in method: MeshPoint::SetAdj() thanks to zuqi.tang@univ-lille.fr)
  -  try to build dmg install mac version
  - | add file script to build meshS from boundary meshL TL if the boundary is the graph of function from  mean plane.
    | see example in examples/3dSurf/buildmeshS.edp

    .. code-block:: freefem
      :linenos:

      meshS Ts=buildmeshSminsurf(TL,1);// minimal surface
      meshS Tsl=buildmeshSLap(TL,1);//  Laplace Surface ..
      meshS Tsl=buildmesh(TL,1,op);// op = 0 Lap and op =1 => minsurf.

  - add sparse block to sparse matrix

    .. code-block:: freefem
      :linenos:

      matrix A = va(Vh,Vh);
      matrix B(A.n*5,A.n*5);
      int i=2;
      B.add(1.+10*i,A,i*ndof,i*ndof);

* Changed

  -  change isoline to do the job for meshS, see example plugin/isoline.edp
  -  change Curve function to be with 3 components to use the isoline data.
  -  change Curvature plugin to compatible with new isoline data for 3 d case.
  -  change some sprintf in snprint to remove warning

* Fixed

  - bug in all P0face, P0edge, P0VF on mesh3,meshS, MeshL  and also discontinous  version (missing  initialisation)
  - bug in  plot function and ffglut with parameter pdf="file.pdf" , because shift in plot named parameter not change in ffglut.
  - genere a bug if zero size element in read MeshL from file.
  - remove mistake when the border is badly defined , remove empty element in buildmeshL function.
  - bug in array quadrature FE.

Version 4.12
------------

* Added

  - | add new finite Element P2pnc3d of Stokes problem like Crouzeix-Raviard in 3d of P2 pylynome
    | see G. Allaire or loic.balaziatchynillama@cea.fr for details
  - | add pdfPLOT from fujiwara@acs.i.kyoto-u.ac.jp (http://www-an.acs.i.kyoto-u.ac.jp/~fujiwara/ff++-programs/)
    | usage: :freefem:`plot( ..., pdf="filename.pdf", svg="filename.svg" );`
  - | add missing code for Discontinous Galerkin in 3d for RHS
    | see `problem-in-3d-discontinuous-galerkin-computation <https://community.freefem.org/t/problem-in-3d-discontinuous-galerkin-computation/2015/6>`__
  - | add in examples/mpi/chamonix.edp : radiative transfer
    | uses new plugin `plugin/mpi/RadiativeTransfer_htool.cpp`, illustrates the use of htool for compression of user defined matrix operator
  - transform a surface meshS in 2d mesh (warning with overlapping, no test) with movemesh:

    .. code-block:: freefem
      :linenos:

      meshS Ths = square3(10,10,[x,y,square(2*x-1)+square(2*y-1)]);
      real[int] gzz;
      mesh Th2 = movemesh(Ths,transfo=[x,y,z],getZ=gzz);//  get flat 2d mesh form meshS

  - New 1d finite element P3 hermite (C1) finite element in plugin `Element_P3`

    .. code-block:: freefem
      :linenos:

      meshL Th=segment(1,[x*L,0,0]); fespace Vh(Th,P3HL);

    see example end of example plugin/testFE-P3
  - missing new 1d finite element P4 in plugin `Element_P4`
  - | plugin `plugin/seq/MatrixMarket.cpp`  to read and save matrix in MatrixMarket and add also a binary form
    | see examples/plugin/MatrixMarket.edp test
  - | add ILU on complex matrix in plugin IncompleteCholesky
    | remark : the IncompleteCholesky is written but not tested
  - add test of functional interface of complex eigen value problem in `examples/eigen/LapEigenValueFuncComplex.edp`

* Changed

  - correct some old code with old version of K.facePermutation() function in plugin/seq/Element_Mixte3d.cpp and plugin/seq/Element_P2bulle3.cpp (not tested)

* Fixed

  - fix in A.RemoveHalf (alway return a new matrix)

Version 4.11
------------

* Added

  - add computation scalar product of R3 example :  ( N'*Tl)
  - add tools to do compution with R3 vector see tutorial/calculus.edp
  - add an example tutorial/tgv-test.edp see see what tgv do on matrix build. 
  - add R3 Th.be(k).N to  get the normal of boundary element (in all mesh type)
  - add R3 Th.be(k)[i].P  to  get the point (R3)  of boundary vertices
  - add R3 Th.be(k).measure to  get the measure of the boundary elment 
  - add projection  function to a mesh , meshL, MeshS or  mesh3 with return a R3 point 
  - see new example dist-projection.edp example in exemples 
  - add dxx, dyy, dzz, dxy,  .. on P2L finite element 
  - add tools to compute solid angle :  let R3 O; a given point, Th3 a mesh3 and ThS a meshS. 
     - solidangle(O,Th3.be(ke)) // triangular face is the boundary face 
     - solidangle(O,Th3[k],nuface) // triangular face is face nuface of tet Th3[k]
     - solidangle(O,ThS[k]) // triangular face is ThS[k]
     - solidangle(O,A,B,C) // triangular face i (A,B,C) 
     - Volume(O,Th3.be(ke)) // O, triangular face is the boundary face 
     - Volume(O,Th3[k],nuface) // O, triangular face is face nuface of tet Th3[k]
     - Volume(O,ThS[k]) // O, triangular face is ThS[k]
     - Volume(O,A,B,C) // (O,A,B,C) tet ..
  - in bem pluging add array of HMatrix    
  -  examples/3d/Connectivite-3d.edp or /3dSurf/Connectivite-S.edp of test. 
  - 3 function mapk, mapkk, mapkk to set a function in fourier space with k parametre

    .. code-block:: freefem
      :linenos:

      R3 K; // le fourier variable allway 3d (sorry)
      int n1=16,n2=8, n3=4; 
      real[int] tab1(nx,tab2(nx*ny),tab3(nx*ny*nz);
      mapk(tab1,K,sqr(K.x));
      mapkk(tab2,ny,K,K.norm2);
      mapkkk(tab3,ny,nz,K,K.norm2);
      //  Remark you can change K by P (current point)
    
  - in SurfaceMesh.ipd fonction to build a Isocaedron and a Sphere from this Isocaedron
  - new finite element on MeshS  this  finite element is the ortogonal of RT0 on surface, or 
    Nelelec Finite Element on triangle with one DoF per mesh edge and where the DoF is the 
    current on  Edge in orientate edge by number of vertices.  
  -  plugin Element_P3pnc for new 2d finite element P3pnc (P3 + 2 bulles)  noncoforming  (continuite of P2 mod)   
      and add 2 examples with this new finite element 
      examples/plugin/cavityNewtowP3pnc.edp examples/plugin/testFE-P3pnc.edp
  - function to set dirichlet Boundary conditon on matrix A (real ou compex) trought  an real[int] 
      (if none zero => set BC ) 
    setBC(A,au1[],-2); and the example 
        examples/3d/Elasticity-simple-support-BC.edp
  
* Changed

  - the beaviour of linear solver UMFPACK, CHOLMOD in case of error , now FreeFEm exit on ExecError like in MUMPS
  - PETSc 3.17.0


* Removed

  -map function  in plugin dfft 

* Fixed

  - pow(int,int) now call int version not complex version..
  - correct the normal the N implicite variable   on meshL case 
  - correct version dump in banner FreeFem++ - version 4.10 (V ...
  - correct  in CPU time on big mesh due to do bad HCode in HashTable.hpp
  - bug in array of finite element on meshhS, meshL (ie.  `fespace Vh(ThS,[P1,P1]);` ) 


Version 4.10
------------

* Added

  - ridgeangle named parameter in ExtractMeshL in msh3 plugin
  - DG formulation in 1d :
    add integral of all border of element : :freefem:`intallBE(ThL)` and unified the notation by adding
    :freefem:`intallBE(ThS)` , :freefem:`intallBE(Th2)`, :freefem:`intallBE(Th3)`
    :freefem:`nuVertex` of now the vertex number of element in :freefem:`intallBE0d` integral
    `BoundaryBE`, `InternalBE` to know if border element (BE) is on true boundary of not.
    update :freefem:`nElementonB` in case on no manifold data (value greater > 2) in meshL, MeshS case ..
    add code to use jump, mean of test functuon on MeshL case. ( not in mesh3 ) to compute RHS.
  - add :freefem:`getcwd()` function in shell plugin to get the current working dir
  - add :freefem:`nuVertex` to get the vextex on element in some int?

* Changed

  - PETSc 3.16.1

* Deprecated

  - SLEPc and SLEPc-complex have been part of PETSc and PETSc-complex for multiple releases and are now deprecated

* Fixed

  - :freefem:`examples/potential.edp` correct problem in times loops and BC
  - :freefem:`tutorial/mortar-DN-4.edp` correct problem of region number in meshL
  - fix problem in Curve mesh and intallBE , vertex number is wrong 
  - portability issue on arm64-apple with `make petsc-slepc`
  - fix assertion failure with `transfer` and `transferMat` with some finite elements


Version 4.9
-----------

* Added

  - add P3 lagrange finite element on meshS and meshS
  - add new plugin :freefem:`meshtool` to add tool to compute the number of connected components of a all kind of mesh
    (mesh,mesh3,meshS,meshL) with 2 kind of connected components ones on interior part of the mesh (default) ans
    secondly on the closure of the mesh (see :freefem:`examples/hpddm/bConnectedComponents.edp` )
    add functions  int[int] In=iminP1K(Th,u) or int[int] Ix=imaxP1K(Th,u)  get the array min/max of value u[i]  
    where i is vertex number on  each element k, so we have  :freefem:`u[Im[k]] = min u[i]/ i in k;`
  - add in plugin `bfstream` to to read binary int (4 bytes) to read fortran file and try to pull tools to share the endiannes
    in progress
  - add gluemesh of array of MeshL and MeshS type
  - interface to :freefem:`PC_MG_GALERKIN_BOTH`
  - Kronecker product of two sparse matrices :freefem:`matrix C = kron(A, B)`
  - add lot of finite element on Mesh3, MeshS, MeshL of Discontinous Galerling Element
    in 3d       : P1dc3d, P2dc3d, P3dc3d, P4dc3d , P0edge3d ,P0edgedc3d ,  P0face3d ,P0facedc3d , P0VF3d ,P0VFdc3d ,
    on Surface  : P1dcS, P2dcS, P3dcS, P4dcS , P0edgeS ,P0edgedcS , P0VFS ,P0VFdcS,
    on Curve   : P1dcL, P2dcL, P3dcL, P4dcL ,  P0VFL ,P0VFdcL
    remark; the associated generic name existe of P1dc, P2dc, P0edge, P0VF and all  dc finite element corresponding to
    no continuity across element.
  - add code of intallfaces to  do Discontinous Galerkin  formulation in 3d (in test FH.)
  - add dist function to a mesh , meshL, MeshS or  mesh3 
  - signeddistfunction to a meshL or  meshS 
  - add buildmesh functon to build a 2d mesh from a meshL (same as buildmesh see examples/3dCurve/border.edp)
 
* Changed

  - Now the order to find MPI in configure is first if you have PETSC then take MPI from PETSc
    otherwise use previous method
  - on MeshL defined with buildmeshL now the default label are 2*k-1  (resp. 2*k)  for the begin (resp. end) of curve
    where k is the order of curve use in buildmeshL. So if you have one curve the  labels are 1  and 2.
    And new  the element label are te region number not the label.
    This element are not really test so be carfull.
  - PETSc 3.15.0


* Fixed

  - bug in Find triangle contening point in 2d (border case),
    :freefem:`int Mesh::DataFindBoundary::Find(R2 PP,R *l,int & outside) const`
    the parameter l not correclty return due to local variable.
  - set CFLAGS=-Wno-implicit-function-declaration to complie with Apple clang version 12.0.0 (clang-1200.0.32.29)
    to remove following error: implicit declaration of function
    correct :freefem:`3dCurve/basicGlue.edp`and add missing test
  - bugs in SLEPc :freefem:`SVDSolve()` with a rectangular :freefem:`Mat`
  - bugs in nElementonB for DG 3d formulation.


Version 4.8
-----------

* Added

  - Bilaplacian example using Morley FE with PETSc, see :freefem:`examples/hpddm/bilaplacian-2d-PETSc.edp`
  - Oseen problem preconditioned by PCD, see :freefem:`examples/hpddm/oseen-2d-PETSc.edp`
  - SLEPc polynomial eigenvalue solver `PEPSolve()`
  - add trivial example to check periodic boundary condition on meshS , meshL  , mesh3
    examples/3d/periodic3.edp	examples/3dSurf/periodicS.edp
    examples/3dCurve/periodicL.edp

* Changed

  - PETSc version 3.14.2
  - Mmg version 5.5.2
  - link of ffglut so change in configure.ac and Makefile.am  LIBS -> FF_LIBS and LIBS become empty
    to remove default libs
  - change number of save plot in ffglut from 10 to 20 for O. Pironneau

* Fixed

  - some memory leaks
  - the periodic boundary condition have wrong before first a sementic level of MeshS and MeshL case.
     the new syntexe is for example:
     meshL Tl=segment(10);   fespace Vl(Tl,P1,periodic=[[1],[2]]);
     meshS Th=square3(10,10,[x*2*pi,y*2*pi]); fespace Vh2(Th,P1,periodic=[[1,x],[3,x],[2,y],[4,y]]);
  - fixed '*' keyboard trick,  to keep  the viewpoint in ffglut or not.


Version 4.7-1
-------------

* Changed

  - change the language definition to use type as a construction function with named arguments for bem plugin
  - PETSc version 3.14.0
  - ARPACK compiled by SLEPc
  - Mmg version 5.5.0
  - -std=c++14 instead of -std=c++11 when possible

* Removed

  - plugins thresholdings, symmetrizeCSR, and fflapack and associed example

* Fixed

  - problem compilation with gfortran-10 of arpack and mumps (add -fallow-argument-mismatch flags)


Version 4.7
-----------

* Added

  - new way to build matrix between 2d Finite element 2d and Curve finite element to do mortar (Thank to Axel ) , see first example `examples/tutorial/mortar-DN-4-v4.5.edp`
  - add :freefem:`Ns` normal vector  in R^3 on meshS (normal of the surface) of current point (to day Ns of [x,y,0] plan  is [0,0,-1])  no be compatible to exterior normal.
  - add :freefem:`Tl` tangent vector in R^3 on meshL (tangent vector of the line/curve) of current point
  - compile ffmaster / ffslave example under windows (thanks to johann@ifado.de)
  - Boolean parameter `spiltpbedge` in `buildmesh` to split in to edge with two boundary vertices
  - interface to PETSc DMPlex, see `examples/hpddm/DMPlex-PETSc.edp`
  - function `MatDestroy`
  - function `MatPtAP` and `transferMat` for parallel interpolation between non-matching grids, see `examples/hpddm/PtAP-2d-PETSc.edp` or `examples/hpddm/diffusion-mg-2d-PETSc.edp`
  - preliminary interface to `SVDSolve` from SLEPc to compute singular value decompositions, see `examples/hpddm/mf-2d-SLEPc.edp` or `examples/hpddm/helmholtz-2d-SLEPc-complex.edp`
  - preliminary interface to `NEPSolve` from SLEPc to solve nonlinear eigenvalue problems, see `examples/hpddm/nonlinear-2d-SLEPc-complex.edp`
  - `transpose` parameter when constructing a `Mat` for defining a matrix-free transposed operation
  - interface to `PetscMemoryGetCurrentUsage`
  - add P2b, RT0, RT1 surface FE (P2bS, RT0S, RT1S))
  - add operator interpolate (2d->3d surface)
  - add operator x = A'\*b; where x, b are array and A 2 dim array (full matrix) and generate an error in case of b'\*A or b'\*A expression
  - function `MatLoad` to load a PETSc `Mat` from disk, see `examples/hpddm/MatLoad-PETSc.edp`
  - possibility to assemble a symmetric `HMatrix<complex>` and to densify a `HMatrix<complex>` into a `Mat<complex>`

* Changed

  - moved Htool to its new GitHub location
  - ScaLAPACK and MUMPS are not compiled by PETSc anymore if there is no Fortran compiler
  - MPICH is compiled by PETSc if no MPI is detected during configure, see https://community.freefem.org/t/feature-request-use-download-mpich-on-ubuntu/407
  - PETSc version 3.13.5
  - force `--with-cudac=0` in `make petsc-slepc`, see https://github.com/FreeFem/FreeFem-sources/issues/141
  - change DSL keyword P1dc3dL->P1dcL and P1dc3dS->P1dcS
  - rename `view`, `hasType`, `changeSchur` to respectively `ObjectView`, `HasType`, and `ChangeSchur`

* Deprecated

  - rename `changeNumbering`, `globalNumbering`, `originalNumbering`, `changeOperator`, `destroyRecycling`, and `attachCoarseOperator` to respectively `ChangeNumbering`, `GlobalNumbering`, `OriginalNumbering`, `ChangeOperator`, `DestroyRecycling`, and `AttachCoarseOperator`
  - `Nt` the normal vector of the current (wrong on meshL) use `Ns` or `Tl`
* Removed 

  - `augmentation` routine from the PETSc plugin
  - `MPIF77` variable

* Fixed

  - lot of mistake in MeshL element add a example o check lot of thing `tutomesh1d.edp`
  - fixed problem of change of mesh when rebuild 2d mesh with buildmesh, .... (Thank to P. Jovilet to points this problem)
  - missing METIS library when using SuiteSparse compiled by PETSc
  - missing `-fno-stack-protector` when building PETSc on Windows, see https://community.freefem.org/t/error-loading-complex-petsc-slepc-library/370
  - fixed ffglut for the plotting of FE array solution
  - fixed  ffglut bug on MacOS Catalina , draw inn only half windows screen (Apple Bug ???)
  - correct P0VF  finite element
  - `abs` function of array


Version 4.6
-----------

* Added

  - new search algorithm for the element containing a point (more safe) in mesh of type :freefem:`mesh3`, :freefem:`meshS`, or :freefem:`meshL`.
  - new function :freefem:`hasType` to know if a PETSc component has been installed, e.g., :freefem:`hasType("PC", "hypre")`
  - eigenvalue problems on linear elements, cf. :freefem:`examples/eigen/LapEigen1DBeltrami.edp` or :freefem:`examples/hpddm/laplace-beltrami-3d-line-SLEPc.edp`
  - `--download-cmake` in PETSc configure if there is no CMake available
  - flags `--with-[slepc|slepccomplex]-include` and `--with-[slepc|slepccomplex]-ldflags` for when SLEPc has been built outside of FreeFEM or PETSc
  - interface to `KSPSetResidualHistory` and `KSPGetIterationNumber`
  - interface to `mpiWaitAll`
  - new function extract, allows to build a curve mesh from a 2d mesh (can extract a labeled boundary, apply a geometric transformation)
  - ffglut can plot a vectorial FE function in surface 3d
  - distributed ParMmg interface, cf. :freefem:`examples/hpddm/distributed-parmmg.edp` or :freefem:`examples/hpddm/laplace-adapt-dist-3d-PETSc.edp`
  - new parallel interpolator on non-matching meshes, cf. :freefem:`examples/hpddm/transfer.edp`
  - ability to solve problems in single precision or with 64 bit integers
  - tool to read data form vtk file only in 3d (cf. plugin iovtk a first example `examples/plugin/iovtk.edp`)
  - tool to read/wrile ply file of meshL, mesh3, MeshS : Polygon File Format / Stanford Triangle Format do  `load "ioply"`
     see :freefem:`examples//3dSurf/operatorsOnMeshS.edp`

* Changed

  - new :freefem:`tgv` values: -10 => zero row, -20 => zero row/column
  - Windows binary now shipped with PETSc/SLEPc
  - BEM examples are now in `examples/mpi`
  - plot border type is now in 3d (border 2d and 3d)
  - PETSc version 3.13.0

* Fixed

  - `--enable-download_package` may now be used to download a single package, e.g., `--enable-download_metis`
  - compilation of PETSc under Windows
  - compilation of plugins when using static libraries
  - correct detection problem in FE type when use a vectorial FE
  - macro concatenation with spaces in arguments
  - correct bug in :freefem:`plugin/seq/Schur-Complement.cpp`
  - correct ambiguity bug in :freefem:`plugin/seq/bfstream.cpp` (reading real or integer)
  - compilation of plugin libff-mmap-semaphore.c under windows


Version 4.5
-----------

Release, binaries packages 
~~~~~~~~~~~~~~~~~~~~~~~~~~

* Since the version 4.5, the FreeFEM binary packages provides with a compiled PETSc library.
* FreeFEM is now interfaced with ParMmg.

New meshes and FEM border 
~~~~~~~~~~~~~~~~~~~~~~~~~

After Surface FEM, Line FEM is possible with a new mesh type :freefem:`meshL`, :freefem:`P0` :freefem:`P1` :freefem:`P2` :freefem:`P1dc` FE, basic FEM, mesh generation.
This new development allows to treat a 1d problem, such as a problem described on a 3d curve.

Abstract about Line FEM in FreeFEM.

* new **meshL** type, refer to the section :ref:`The type meshL in 3 dimension <meshStype>`  
  
  - new type of surface mesh: :freefem:`meshL`
  - the functionalities on the :freefem:`meshL` type, it is necessary to load the plugin ”msh3”. 
  - generator of meshL :freefem:`segment`, define multi :freefem:`border` and :freefem:`buildmesh` function.
  - basic transformation are avalaible: :freefem:`movemesh`, :freefem:`trunc`, :freefem:`extract`, :freefem:`checkmesh`, :freefem:`change`, :freefem:`AddLayers`, glue of :freefem:`meshL`.
  
    It is possible to build the underlying :freefem:`meshL` from a :freefem:`meshS` with the function :freefem:`buildBdMesh`: :freefem:`ThS=buildBdMesh(ThS)` builds the boundary domain associated to the :freefem:`meshS` ThS and extract it by the command :freefem:`meshL` ThL=ThS. :freefem:`Gamma`. 
  

* new finite element space with curve finite element type
 
 - :freefem:`FESpace` :freefem:`P0` :freefem:`P1`, :freefem:`P2`, :freefem:`P1dc` Lagrange finite elements and possible to add a custumed finite element with the classical method (like a plugin).

* as in the standard 2d, 3d, surface 3d case, the variational problem associated to surface PDE can be defined by using the keywords

  - :freefem:`problem` 
  - :freefem:`varf` to access to matrix and RHS vector
  - available operators are :freefem:`int1d`, :freefem:`on` and the operator :freefem:`int0d` to define a Neumann boundary condition 


* visualisation tools 

  - plot with :freefem:`plot` of ffglut, :freefem:`medit` meshes meshL and solutions
  - 2d or 3d view, with in 3d the option to visualize the elememt Normals at element (touch 'T') and the deformed domain according to it (touch '2').
  - loading, saving of meshes and solution at FreeFEM's format
    
    + ".mesh"  mesh format file of Medit (P. Frey LJLL) 
    + ".msh" for mesh and ".sol" data solution at freefem format
    + ".msh" data file of Gmsh (Mesh generator) (load  "gmsh")
    + vtk format for meshes and solutions (load "iovtk" and use the ".vtu" extension)


Boundary Element Method
~~~~~~~~~~~~~~~~~~~~~~~

Allows to define and solve a 2d/3d BEM formulation and rebuild the associated potential.
The document is in construction.
