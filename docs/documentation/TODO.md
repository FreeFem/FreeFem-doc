<!--- THIS FILE IS AUTOMATICALY GENERATED --->
<!--- DO NOT EDIT --->

# TODO

## Home

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Notations

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Mesh generation

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:99%">
	</div>
	<span class="progress-label">99</span>
</div>

- [ ] line 2597
	``` mmg3d-v4.0```

## Finite element

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Visualization

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Parallelization

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:98%">
	</div>
	<span class="progress-label">98</span>
</div>

- [ ] line 54
	``` script bug```
- [ ] line 132
	``` script bug```
- [ ] line 154
	``` mpiReduceScatter is commented in parallelempi.cpp```
- [ ] line 156
	```See the `:::freefem examples++-mpi/essai.edp`  to test of all this functionality and thank you to Guy-Antoine Atenekeng Kahou for his help to code this interface.```
- [ ] line 237
	``` script bug```
- [ ] line 241
	``` check this part```

## Chapter 11

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:94%">
	</div>
	<span class="progress-label">94</span>
</div>

- [ ] line 28
	```If the libraries are not loaded, the default sparse solver will be loaded (default sparse solver is UMFPACK). The table 11.1 \ref{lib.sparse.solver}  gives this new value for the different libraries.```
- [ ] line 88
	```We also add functions (see Table 11.2 \ref{func.sparse.solver} ) with no parameter to change the default sparse solver in the .edp file. To use these functions, we need to load the library corresponding to the solver. An example of using different parallel sparse solvers for the same problem is given in testdirectsolvers.edp (directory example$++-$mpi ).```
- [ ] line 238
	``````
- [ ] line 252
	```MUltifrontal Massively Parallel Solver (MUMPS) is a free library \cite{mumpspubl1,mumpspubl2,mumpspubl3} . This package solves linear system of the form $A \: x = b$ where $A$ is a square sparse matrix with a direct method. The square matrix considered in MUMPS can be either unsymmetric, symmetric positive definite or general symmetric.```
- [ ] line 253
	```The method implemented in MUMPS is a direct method based on a multifrontal approach \cite{mumpspubl1} . It constructs a direct factorization $A \:= \: L\:U$, $A\: = \: L^t \: D \: L$ depending of the symmetry of the matrix $A$. MUMPS uses the following libraries : BLAS\cite{blas1,blas2} , BLACS and ScaLAPACK\cite{scalapackuserguide} .```
- [ ] line 263
	```MUMPS is written in Fortran 90. The parallel version is constructed using MPI \cite{mpi}  for message passing and BLAS \cite{blas1,blas2} , BLACS and ScaLAPACK\cite{scalapackuserguide} . Therefore, a fortran compiler is needed, and MPI, BLAS, BLACS and ScaLAPACK. An installation procedure to obtain this package is given in the file README\_COMPILE in the directory src/solver  of FreeFem++.```
- [ ] line 267
	```The MUMPS interface for FreeFem++ is given in file MUMPS\_freefem.cpp (directory src/solver/ ).```
- [ ] line 269
	```A description to obtain this library is given in the file README\_COMPILE in the directory src/solver  of FreeFem++. We recall here the procedure. Go to the directory src/solver in FreeFem++ package. Edit the file makefile-sparsesolver.inc to yours system: comment Section 1, comment line corresponding to libraries BLAS, BLACS, ScaLAPACK, Metis, scotch in Section 2 and comment in Section 3  the paragraph corresponding to MUMPS solver. And then type `:::bash make mumps` in a terminal window.```
- [ ] line 275
	```There are four input parameters in MUMPS (see \cite{mumpsuserguide} ). Two integers SYM and PAR, a vector of integer of size 40 INCTL and a vector of real of size 15 CNTL. The first parameter gives the type of the matrix: 0 for unsymmetric matrix, 1 for symmetric positive matrix and 2 for general symmetric. The second parameter defined if the host processor work during the factorization and solves steps : PAR=1 host processor working and PAR=0 host processor not working.```
- [ ] line 276
	```The parameter INCTL and CNTL is the control parameter of MUMPS. The vectors ICNTL and CNTL in MUMPS becomes with index 1 like vector in fortran. A short description of all parameters of ICNTL and CNTL is given in ffmumps\_fileparam.txt . For more details see the users' guide \cite{mumpsuserguide} .```
- [ ] line 290
	```	where $P$ is the permutation matrix, $Q_c$ is the column permutation, $D_r$ and $D_c$ are diagonal matrix for respectively row and column scaling. The ordering strategy to obtain $P$ is controlled by parameter ICNTL(7). The permutation of zero free diagonal $Q_c$ is controlled by parameter ICNTL(6). The row and column scaling is controlled by parameter ICNTL(18). These option are connected and also strongly related with ICNTL(12) (see documentation of mumps for more details \cite{mumpsuserguide} ). The parameters `:::freefem permr`, `:::freefem scaler`, and `:::freefem scalec` in FreeFem++ allow to give permutation matrix($P$), row scaling ($D_r$) and column scaling ($D_c$) of the user respectively.```
- [ ] line 294
	```To call MUMPS in FreeFem++, we need to load the dynamic library MUMPS\_freefem.dylib (MacOSX), MUMPS\_freefem.so (Unix) or MUMPS\_freefem.dll (Windows) .```
- [ ] line 373
	```A simple example of calling MUMPS in FreeFem++ with this two methods is given in the file testsolver_MUMPS.edp  in the directory examples++-mpi.```
- [ ] line 377
	```The package SuperLU_DIST \cite{slu2,slu1}  solves linear systems using LU factorization. It is a free scientific library under BSD license. The web site of this project is [http://crd.lbl.gov/~xiaoye/SuperLU](http://crd.lbl.gov/~xiaoye/SuperLU). This library provides functions to handle square or rectangular matrix in real and complex arithmetics. The method implemented in SuperLU_DIST is a supernodal method \cite{slu1} . New release of this package includes a parallel symbolic factorization \cite{slu2} . This scientific library is written in C and MPI for communications.```
- [ ] line 381
	```To use SuperLU_DIST in FreeFem++, you have to install SuperLU_DIST package. We need MPI and ParMetis library to do this compilation. An installation procedure to obtain this package is given in the file README\_COMPILE in the directory src/solver/  of the FreeFem++ package.```
- [ ] line 385
	```The FreeFem++ interface to SuperLU_DIST for real (resp. complex) arithmetics is given in file real_SuperLU_DIST_FreeFem.cpp (resp. complex_SuperLU_DIST_FreeFem.cpp). These files are in the directory src/solver/. These interfaces are compatible with the release 3.2.1 of SuperLU_DIST. To use SuperLU_DIST in FreeFem++, we need libraries corresponding to these interfaces. A description to obtain these libraries is given in the file README_COMPILE in the directory src/solver  of FreeFem++. We recall here the procedure. Go to the directory src/solver in FreeFem++ package. Edit the file makefile-sparsesolver.inc in your system : comment Section 1, comment line corresponding to libraries BLAS, Metis, ParMetis in Section 2 and comment in Section 3 the paragraph corresponding to SuperLU_DIST solver. And just type `:::bash make rsludist` (resp. `:::bash make csludist`) in the terminal to obtain the dynamic library of interface for real (resp. complex) arithmetics.```
- [ ] line 399
	```The option arguments of SuperLU_DIST are described in the section Users-callable routine of \cite{sluuserguide} . The parameter Fact and TRANS are specified in FreeFem++ interfaces to SuperLU_DIST during the different steps. For this reason, the value given by the user for this option is not considered.```
- [ ] line 408
	```The other parameters for LU factorization are ParSymFact and ReplaceTinyPivot. The parallel symbolic factorization works only on a power of two processes and need the ParMetis ordering \cite{parmetis} . The default option argument of SuperLU_DIST are given in the file ffsuperlu_dist_fileparam.txt.```
- [ ] line 435
	```This value correspond to the parameter in the file ffsuperlu_dist_fileparam.txt.  If one parameter is not specified by the user, we take the default value of SuperLU_DIST.```
- [ ] line 463
	```A simple example of calling SuperLU_DIST in FreeFem++ with this two methods is given in the file testsolver_superLU_DIST.edp in the directory examples++-mpi. ```
- [ ] line 467
	```Pastix (Parallel Sparse matrix package) is a free scientific library under CECILL-C license. This package solves sparse linear system with a direct and block ILU(k) iterative methods. This solver can be applied to a real or complex matrix with a symmetric pattern \cite{pastix} .```
- [ ] line 471
	```To used Pastix in FreeFem++, you have to install pastix package in first. To compile this package, we need a fortran 90 compiler, scotch \cite{scotch}  or Metis \cite{metis}  ordering library and MPI. An installation procedure to obtain this package is given in the file .src/solver/ README_COMPILE in the section pastix of the FreeFem++ package.```
- [ ] line 481
	```The input `:::freefem matrix` parameter of FreeFem++ depend on pastix interface. `:::freefem matrix = assembled` for non distributed matrix. It is the same parameter for SuperLU_DIST. There are four parameters in Pastix : iparm, dparm, perm and invp. These parameters are respectively the integer parameters (vector of size 64), real parameters (vector of size 64), permutation matrix and inverse permutation matrix respectively. iparm and dparm vectors are described in \cite{pastixrefcard} .```
- [ ] line 487
	``````
- [ ] line 517
	```A simple example of calling pastix in FreeFem++ with this two methods is given in the file testsolver_pastix.edp in the directory examples++-mpi. ```
- [ ] line 519
	```In Table 11.3 \ref{recap.direct.solveur} , we recall the different matrix considering in the different direct solvers.```
- [ ] line 574
	```Concerning __iterative solvers__, we have chosen _pARMS_ \cite{spARMS} , _HIPS_ \cite{HIPS}  and _Hypre_ \cite{HYPRE} .```
- [ ] line 576
	```So, _pARMS_ implements algebraic domain decomposition preconditioner type such as additive Schwartz \cite{CAI89}  and interface method \cite{LISAAD} ; while HIPS implement hierarchical incomplete factorization \cite{A:LaBRI::sisc06}  and finally HYPRE implements multilevel preconditioner are AMG(Algebraic MultiGrid) \cite{A:AMG::sisc00}  and parallel approximated inverse \cite{A:PASA::sisc00} .```
- [ ] line 581
	```METIS \cite{metis}  or Scotch \cite{scotch} .```
- [ ] line 590
	```_pARMS_ (parallel Algebraic Multilevel Solver) is a software developed by Youssef Saad and al at University of Minnesota \cite{spARMS} .```
- [ ] line 592
	```It consists of variants of GMRES like FGMRES (Flexible GMRES) , DGMRES (Deflated GMRES) \cite{SAAD03}  and BICGSTAB. pARMS also implements parallel preconditioner like RAS (Restricted Additive Schwarz)\cite{CAI89}  and Schur Complement type preconditioner \cite{LISAAD} .```
- [ ] line 594
	```All these parallel preconditioners are based on the principle of domain decomposition. Thus, the matrix $A$ is partitioned into sub matrices $A_i$($i=1,...,p$) where p represents the number of partitions one needs. The union of $A_i$ forms the original matrix. The solution of the overall system is obtained by solving the local systems on $A_i$ (see \cite{Smith96} ).```
- [ ] line 600
	```To install pARMS, you must first download the pARMS package at \cite{spARMS} . Once the download is complete, you must unpack package pARMS and follow the installation procedure described in file README to create the library __libparms.a__.```
- [ ] line 603
	```Before calling pARMS solver inside FreeFem++, you must compile file parms\_FreeFem.cpp to create a dynamic library parms\_FreeFem.so. To do this, move to the directory src/solver of FreeFem++, edit the file makefile\-parms.inc to specify the following variables: ```
- [ ] line 623
	```This example comes from user guide of FreeFem++ \cite{ufreefem}  at page 12.```
- [ ] line 643
	```In line 1 of example \ref{exm:first}  we load in memory the pARMS dynamic library```
- [ ] line 654
	```* Tolerance for convergence=$1e-08$.(see book of Saad \cite{SAAD03}  to```
- [ ] line 656
	```* preconditionner=Restricted Additif Schwarz \cite{CAI89} ,```
- [ ] line 667
	```Lets us consider Navier Stokes example \ref{exm:segond}  . In this example we```
- [ ] line 697
	```We need two vectors to specify the parameters of the linear solver. In line 1 of example \ref{exm:segond}  we have declared these vectors(`:::freefem int[int] iparm(16); real[int] dparm(6);`). In line 3 we have initialized these vectors by negative values. We do this because all parameters values in pARMS are positive and if you do not change the negative values of one entry of this vector, the default value will be set. In tables (table 11.4 \ref{lpparm}  and 11.5 \ref{pardoubleparm} ), we have the meaning of differents entries of these vectors.```
- [ ] line 703
	```			\ref{exm:segond} </th>```
- [ ] line 713
	```	        <td>Krylov subspace methods.<br>Differents valuesfor this parameters are specify on table \ref{kryparms} </td>```
- [ ] line 717
	```	        <td>Preconditionner.<br>Differentspreconditionners for this parameters are specify on table \ref{precond} </td>```
- [ ] line 785
	```			<th colspan="2">Table 11.5 : Significations of dparams corresponding variables for example \ref{exm:segond} </th>```
- [ ] line 916
	```We run example 11.4 \ref{exm:segond}  on cluster paradent of Grid5000 and report results in table 11.8 \ref{parmResult} .```
- [ ] line 922
	```			\ref{exm:segond} 11.4 </th>```

## Chapter 11 part2

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:98%">
	</div>
	<span class="progress-label">98</span>
</div>

- [ ] line 3
	```__HIPS__ (_Hierarchical Iterative Parallel Solver_) is a scientific library that provides an efficient parallel iterative solver for very large sparse linear systems. HIPS is available as free software under the CeCILL-C licence. The interface that we realized is compatible with release __1.2 beta.rc4__  of HIPS.```
- [ ] line 5
	```HIPS implements two solver classes which are the iteratives class (GMRES, PCG) and the Direct class. Concerning preconditionners, HIPS implements a type of multilevel ILU. For further informations on those preconditionners see \cite{A:LaBRI::sisc06,A:LaBRI::HRR07} .```
- [ ] line 9
	```To install HIPS, first download the HIPS package at \cite{HIPS} , unpack it and go to the HIPS source directory. The installation of HIPS is machine dependence. For example, to install HIPS on a linux cluster copy the file __Makefile\_Inc\_Examples/makefile.inc.gnu__  on the root directory of HIPS with the name __makefile.inc__. After this, edit __makefile.inc__ to set values of different variables and type `:::bash make all`.```
- [ ] line 13
	```Before calling the HIPS solver inside FreeFem++, you must compile file hips\_FreeFem.cpp to create dynamic library hips\_FreeFem.so. To do this, move to the directory src/solver of FreeFem++ and edit the file makefile.inc to specify the following variables: ```
- [ ] line 17
	```HIPS_INCLUDE : -I$(HIPS_DIR)/SRC/INCLUDE : Directory for HIPS headers ```
- [ ] line 18
	```LIB_DIR : -L$(HIPS_DIR)/LIB : Librairies directory ```
- [ ] line 34
	```Let us consider the 3D Laplacian example inside FreeFem++ package where after discretization we want to solve the linear equation with Hips. Example \ref{hips:laplacian}  is Laplacian3D using Hips as linear solver. We first load Hips solver at line 2. From line 4 to 15 we specify the parameters for the Hips solver and in line 46 of example \ref{hips:laplacian}  we set these parameters in the linear solver.```
- [ ] line 36
	```In Table 1.10 \ref{hipslabel}  results of running example 1.1 5\ref{hips:laplacian}  on Cluster Paradent of Grid5000 are reported. We can see in this running example the efficiency of parallelism.```
- [ ] line 101
	```			<th colspan="3">Table 11.10 : Legend of table 11.10 are give in table 11.9 \ref{legtableparm} </th>```
- [ ] line 253
	```__HYPRE__ (High Level Preconditioner) is a suite of parallel preconditioner developed at Lawrence Livermore National Lab \cite{HYPRE}  .```
- [ ] line 260
	```The _sparse approximate inverse_ approximates the inverse of a matrix $A$ by a sparse matrix $M$. A technical idea to construct matrix $M$ is to minimize the Frobenuis norm of the residual matrix $I-MA$. For more details on this preconditioner technics see \cite{chow} .```
- [ ] line 266
	```To install HYPRE, first download the HYPRE package at \cite{HYPRE} , unpack it and go to the HYPRE/src source directory and do `:::bash ./configure` to configure Hypre. After this just type `:::bash make all` to create __libHYPRE.a__.```
- [ ] line 290
	```Let us consider again the 3D Laplacian example inside FreeFem++ package where after discretization we want to solve the linear equation with Hypre. Example \ref{hypre:laplacian}  is the Laplacian3D using Hypre as linear solver. Example \ref{hypre:laplacian}  is the same as \ref{hips:laplacian} , so we just show here the lines where we set some Hypre parameters.```
- [ ] line 294
	```It should be noted that the meaning of the entries of these vectors is different from those of Hips. In the case of HYPRE, the meaning of differents entries of vectors `:::freefem iparm` and `:::freefem dparm` are given in tables \ref{communipramsHypre}  to \ref{AMGHypre} .```
- [ ] line 296
	```In Table \ref{hyprelabel}  the results of running example \ref{hypre:laplacian}  on Cluster Paradent of Grid5000 are reported. We can see in this running example the efficiency of parallelism, in particular when AMG are use as preconditioner.```
- [ ] line 596
	```\ref{exm:segond}  }```

## Plugins

Progression:
<div class="progress progress-100plus">
	<div class="progress-bar" style="width:100%">
	</div>
	<span class="progress-label">100</span>
</div>


## Developers

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:98%">
	</div>
	<span class="progress-label">98</span>
</div>

- [ ] line 596
	``````
- [ ] line 663
	```It is agood idea to first try the example `load.edp` in directory `example++-load` .```
- [ ] line 671
	```Now, assume that you are in a shell window (a `cygwin` window under Windows) in the directory `example++-load` .```
- [ ] line 807
	```This will add FFT to __`FreeFem++`__, taken from \url{http://www.fftw.org/}. To download and install under `download/include` just go in `download/fftw` and try `make`.```
- [ ] line 850
	```To test, try `dfft.edp`. ```
- [ ] line 854
	```The associed edp file is `examples++-load/convect_dervieux.edp`. ```
- [ ] line 856
	```See `mat_dervieux.cpp`. ```
- [ ] line 860
	```First read the [Adding a new finite element section](#adding-a-new-finite-element), we add two new finite elements examples in the directory `examples++-load`. ```
- [ ] line 878
	```See `BernardiRaugel.cpp`. ```
- [ ] line 927
	```A real example using this finite element, just a small modification of the `NSP2P1.edp`  examples, just the begenning is change to```
- [ ] line 942
	``` PR needed: BernaRdiRaugel.cpp```
- [ ] line 945
	```See the example `bilapMorley.edp` .```
- [ ] line 950
	```Warning the sparse solver interface as been completely rewritten in version 3.2, so the section is obsolete, the example in are correct/ ```

