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
	<div class="progress-bar" style="width:99%">
	</div>
	<span class="progress-label">99</span>
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

## Chapter 11 part2

Progression:
<div class="progress progress-80plus">
	<div class="progress-bar" style="width:97%">
	</div>
	<span class="progress-label">97</span>
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
- [ ] line 537
	```			<th colspan="3">Table 11.18 : Convergence and time for solving linear system from example 11.4 \ref{exm:segond}  </th>```
- [ ] line 698
	````:::freefem mpiReduce(Data a,Data b,processor(int rk, mpiComm cc),MPI_Op op)`,  ERREUR DE FRAPPE QQ PART ? Reduces values `:::freefem Data a` on all processes to a single value `:::freefem Data b` on process of rank `:::freefem rk` and communicator `:::freefem cc`.```
- [ ] line 711
	``````
- [ ] line 848
	```	* `:::freefem D` is the diagonal of the local partition of unity (see below \S~\ref{sub:linear} 11.5.2  for more details)```
- [ ] line 860
	```In the above line, the first option selects the one-level preconditioner `:::freefem ras` (possible choices are `:::freefem ras`, `:::freefem oras`, `:::freefem soras`, `:::freefem asm`, `:::freefem osm` or `:::freefem none`), the second option selects the correction formula for the second level here `:::freefem balanced` (possible options are `:::freefem deflated`, `:::freefem additive` or `:::freefem balanced`), the third option selects right preconditioning, the fourth one is verbosity level of HPDDM (different from the one of FreeFem++), the fifth one prints all possible options of HPPDM and the last one specifies the number of coarse degrees of freedom per subdomain of the GENEO coarse space. All other options of [https://github.com/hpddm/hpddm/blob/master/doc/cheatsheet.pdf](cheatsheet of the HPDDM) \cite{Jolivet:2014:HPD}  library can be selected via the FreeFem++ function `:::freefem set`.```
- [ ] line 1026
	```As an example, consider the scalar product of two distributed vectors ${\mathbf U}, {\mathbf V} \in \mathbb{R}^{n}$. Using the partition of unity~\eqref{eq:hpddm:14} , we have:```
- [ ] line 1046
	```The matrix vector product is more involved and details are given in the SIAM book  [https://www.ljll.math.upmc.fr/nataf/OT144DoleanJolivetNataf_full.pdf](An Introduction to Domain Decomposition Methods: algorithms, theory and parallel implementation) \cite{Dolean:2015:IDD}  and even more details are given in [http://jolivet.perso.enseeiht.fr/thesis.pdf](P. Jolivet's PhD manuscrit).```

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

