# Parallel sparse solvers

Parallel sparse solvers use several processors to solve linear systems of equation. Like sequential, parallel linear solvers can be direct or iterative. In __`FreeFem++`__ both are available.

## Using parallel sparse solvers in __`FreeFem++`__

We recall that the `:::freefem solver` parameters are defined in the following commands: `:::freefem solve`, `:::freefem problem`, `:::freefem set` (setting parameter of a matrix) and in the construction of the matrix corresponding to a bilinear form. In these commands, the parameter `:::freefem solver` must be set to `:::freefem sparsesolver` for parallel sparse solver. We have added specify parameters to these command lines for parallel sparse solvers. These are

* `:::freefem lparams` : vector of integer parameters (`l` is for the `C++` type `:::cpp long`)
* `:::freefem dparams` : vector of real parameters
* `:::freefem sparams` : string parameters
* `:::freefem datafilename` : name of the file which contains solver parameters

The following four parameters are only for direct solvers and are vectors. These parameters allow the user to preprocess the matrix (see the section on [sparse direct solver](#sparse-direct-solver) for more information).

* `:::freefem permr` : row permutation (integer vector)
* `:::freefem permc` : column permutation or inverse row permutation (integer vector)
* `:::freefem scaler` : row scaling (real vector)
* `:::freefem scalec` : column scaling (real vector)

There are two possibilities to control solver parameters. The first method defines parameters with `:::freefem lparams`, `:::freefem dparams` and `:::freefem sparams` in `.edp` file.

The second one reads the solver parameters from a data file. The name of this file is specified by `:::freefem datafilename`. If `:::freefem lparams`, `:::freefem dparams`, `:::freefem sparams` or `:::freefem datafilename` is not provided by the user, the solver's default values are used.

To use parallel solver in __`FreeFem++`__, we need to load the dynamic library corresponding to this solver. For example to use [MUMPS](http://mumps.enseeiht.fr/) solver as parallel solver in FreeFem, write in the `.edp` file `:::freefem load "MUMPS_FreeFem"`.

If the libraries are not loaded, the default sparse solver will be loaded (default sparse solver is `:::freefem UMFPACK`). The [table 1](#Tab1) gives this new value for the different libraries.

<table>
	<thead>
		<tr>
			<th colspan="3"><a name="Tab1">Table 1</a>: Default sparse solver for real and complex arithmetics when we load a parallel sparse solver library</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td rowspan="2" style="vertical-align: middle !important; font-weight: bold">Libraries</td>
			<td colspan="2" align="center" style="font-weight: bold">Default sparse solver</td>
		</tr>
		<tr>
			<td align="center" style="font-weight: bold">real</td>
			<td align="center" style="font-weight: bold">complex</td>
		</tr>
		<tr>
			<td>MUMPS_FreeFem</td>
			<td align="center">mumps</td>
			<td align="center">mumps</td>
		</tr>
		<tr>
			<td>real_SuperLU_DIST_FreeFem</td>
			<td align="center">SuperLU_DIST</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>complex_SuperLU_DIST_FreeFem</td>
			<td align="center">previous solver</td>
			<td align="center">SuperLU_DIST</td>
		</tr>
		<tr>
			<td>real_pastix_FreeFem</td>
			<td align="center">pastix</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>complex_pastix_FreeFem</td>
			<td align="center">previous solver</td>
			<td align="center">pastix</td>
		</tr>
		<tr>
			<td>hips_FreeFem</td>
			<td align="center">hips</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>hypre_FreeFem</td>
			<td align="center">hypre</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>parms_FreeFem</td>
			<td align="center">parms</td>
			<td align="center">previous solver</td>
		</tr>
	</tbody>
</table>

We also add functions (see [Table 2](#Tab2)) with no parameter to change the default sparse solver in the `.edp` file. To use these functions, we need to load the library corresponding to the solver. An example of using different parallel sparse solvers for the same problem is given in [Direct solvers example](../examples/#direct-solvers).

<table>
	<thead>
		<tr>
			<th colspan="3"><a name="Tab2">Table 2</a>: Functions that allow to change the default sparse solver for real and complex arithmetics and the result of these functions</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td rowspan="2" style="vertical-align: bottom !important; font-weight: bold">Function</td>
			<td colspan="2" align="center" style="font-weight: bold">default sparse solver</td>
		</tr>
		<tr>
			<td align="center" style="font-weight: bold">real</td>
			<td align="center" style="font-weight: bold">complex</td>
		</tr>
		<tr>
			<td>defaulttoMUMPS()</td>
			<td align="center">mumps</td>
			<td align="center">mumps</td>
		</tr>
		<tr>
			<td>realdefaulttoSuperLUdist()</td>
			<td align="center">SuperLU_DIST</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>complexdefaulttoSuperLUdist()</td>
			<td align="center">previous solver</td>
			<td align="center">SuperLU_DIST</td>
		</tr>
		<tr>
			<td>realdefaultopastix()</td>
			<td align="center">pastix</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>complexdefaulttopastix()</td>
			<td align="center">previous solver</td>
			<td align="center">pastix</td>
		</tr>
		<tr>
			<td>defaulttohips()</td>
			<td align="center">hips</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>defaulttohypre()</td>
			<td align="center">hypre</td>
			<td align="center">previous solver</td>
		</tr>
		<tr>
			<td>defaulttoparms()</td>
			<td align="center">parms</td>
			<td align="center">previous solver</td>
		</tr>
	</tbody>
</table>

!!!question "Test direct solvers"
	```freefem
	load "MUMPS_FreeFem"
	//default solver: real-> MUMPS, complex -> MUMPS
	load "real_SuperLU_DIST_FreeFem"
	//default solver: real-> SuperLU_DIST, complex -> MUMPS
	load "real_pastix_FreeFem"
	//default solver: real-> pastix, complex -> MUMPS

	// Solving with pastix
	{
		matrix A =
			[[1, 2, 2, 1, 1],
			[ 2, 12, 0, 10 , 10],
			[ 2, 0, 1, 0, 2],
			[ 1, 10, 0, 22, 0.],
			[ 1, 10, 2, 0., 22]];

		real[int] xx = [1, 32, 45, 7, 2], x(5), b(5), di(5);
		b = A*xx;
		cout << "b = " << b << endl;
		cout << "xx = " << xx << endl;

		set(A, solver=sparsesolver, datafilename="ffpastix_iparm_dparm.txt");
		cout << "solve" << endl;
		x = A^-1*b;
		cout << "b = " << b << endl;
		cout << "x = " << endl;
		cout << x << endl;
		di = xx - x;
		if (mpirank == 0){
			cout << "x-xx = " << endl;
			cout << "Linf = " << di.linfty << ", L2 = " << di.l2 << endl;
		}
	}

	// Solving with SuperLU_DIST
	realdefaulttoSuperLUdist();
	//default solver: real-> SuperLU_DIST, complex -> MUMPS
	{
		matrix A =
			[[1, 2, 2, 1, 1],
			[ 2, 12, 0, 10 , 10],
			[ 2, 0, 1, 0, 2],
			[ 1, 10, 0, 22, 0.],
			[ 1, 10, 2, 0., 22]];

		real[int] xx = [1, 32, 45, 7, 2], x(5), b(5), di(5);
		b = A*xx;
		cout << "b = " << b << endl;
		cout << "xx = " << xx << endl;

		set(A, solver=sparsesolver, datafilename="ffsuperlu_dist_fileparam.txt");
		cout << "solve" << endl;
		x = A^-1*b;
		cout << "b = " << b << endl;
		cout << "x = " << endl;
		cout << x << endl;
		di = xx - x;
		if (mpirank == 0){
			cout << "x-xx = " << endl;
			cout << "Linf = " << di.linfty << ", L2 = " << di.l2 << endl;
		}
	}

	// Solving with MUMPS
	defaulttoMUMPS();
	//default solver: real-> MUMPS, complex -> MUMPS
	{
		matrix A =
			[[1, 2, 2, 1, 1],
			[ 2, 12, 0, 10 , 10],
			[ 2, 0, 1, 0, 2],
			[ 1, 10, 0, 22, 0.],
			[ 1, 10, 2, 0., 22]];

		real[int] xx = [1, 32, 45, 7, 2], x(5), b(5), di(5);
		b = A*xx;
		cout << "b = " << b << endl;
		cout << "xx = " << xx << endl;

		set(A, solver=sparsesolver, datafilename="ffmumps_fileparam.txt");
		cout << "solving solution" << endl;
		x = A^-1*b;
		cout << "b = " << b << endl;
		cout << "x = " << endl;
		cout << x << endl;
		di = xx - x;
		if (mpirank == 0){
			cout << "x-xx = " << endl;
			cout << "Linf = " << di.linfty << ", L2 " << di.l2 << endl;
		}
	}
	```

## Sparse direct solver

In this section, we present the sparse direct solvers interfaced with __`FreeFem++`__.

### MUMPS solver

MUltifrontal Massively Parallel Solver ([MUMPS](http://mumps.enseeiht.fr/)) is an open-source library.

This package solves linear system of the form $A \: x = b$ where $A$ is a square sparse matrix with a direct method. The square matrix considered in MUMPS can be either unsymmetric, symmetric positive definite or general symmetric.

The method implemented in MUMPS is a direct method based on a multifrontal approach. It constructs a direct factorization $A \:= \: L\:U$, $A\: = \: L^t \: D \: L$ depending of the symmetry of the matrix $A$.

MUMPS uses the following libraries : [BLAS](http://www.netlib.org/blas/), [BLACS](http://www.netlib.org/blacs/) and [ScaLAPACK](http://www.netlib.org/scalapack/).

!!!warning
	MUMPS does not solve linear system with a rectangular matrix.

__MUMPS parameters:__

There are four input parameters in [MUMPS](http://mumps.enseeiht.fr/index.php?page=doc). Two integers `:::cpp SYM` and `:::cpp PAR`, a vector of integer of size 40 `:::cpp INCTL` and a vector of real of size 15 `:::cpp CNTL`.

The first parameter gives the type of the matrix: 0 for unsymmetric matrix, 1 for symmetric positive matrix and 2 for general symmetric.

The second parameter defined if the host processor work during the factorization and solves steps : `:::cpp PAR=1` host processor working and `:::cpp PAR=0` host processor not working.

The parameter `:::cpp INCTL` and `:::cpp CNTL` is the control parameter of MUMPS. The vectors `:::cpp ICNTL` and `:::cpp CNTL` in MUMPS becomes with index 1 like vector in fortran. For more details see the [MUMPS users' guide](http://mumps.enseeiht.fr/index.php?page=doc).

We describe now some elements of the main parameters of `:::cpp ICNTL` for MUMPS.

* __Input matrix parameter__
	The input matrix is controlled by parameters `ICNTL(5)` and `ICNTL(18)`. The matrix format (resp. matrix pattern and matrix entries) are controlled by `INCTL(5)` (resp. `INCTL(18)`).

	The different values of `ICNTL(5)` are 0 for assembled format and 1 for element format. In the current release of __`FreeFem++`__, we consider that FE matrix or matrix is storage in assembled format. Therefore, `INCTL(5)` is treated as 0 value.

	The main option for `ICNTL(18)`: `INCLTL(18)=0` centrally on the host processor, `ICNTL(18)=3` distributed the input matrix pattern and the entries (recommended option for distributed matrix by developer of MUMPS). For other values of `ICNTL(18)` see the [MUMPS users' guide](http://mumps.enseeiht.fr/index.php?page=doc). These values can be used also in __`FreeFem++`__.

	The default option implemented in __`FreeFem++`__ are `ICNTL(5)=0` and `ICNTL(18)=0`.

* __Preprocessing parameter__
	The preprocessed matrix $A_{p}$ that will be effectively factored is defined by
	$$
	A_{p} = P \: D_r \: A \: Q_c \ D_c P^t
	$$
	where $P$ is the permutation matrix, $Q_c$ is the column permutation, $D_r$ and $D_c$ are diagonal matrix for respectively row and column scaling.

	The ordering strategy to obtain $P$ is controlled by parameter `ICNTL(7)`. The permutation of zero free diagonal $Q_c$ is controlled by parameter `ICNTL(6)`. The row and column scaling is controlled by parameter `ICNTL(18)`. These option are connected and also strongly related with `ICNTL(12)` (see the [MUMPS users' guide](http://mumps.enseeiht.fr/index.php?page=doc) for more details).

	The parameters `:::freefem permr`, `:::freefem scaler`, and `:::freefem scalec` in __`FreeFem++`__ allow to give permutation matrix($P$), row scaling ($D_r$) and column scaling ($D_c$) of the user respectively.

__Calling MUMPS in `FreeFem++`__

To call MUMPS in __`FreeFem++`__, we need to load the dynamic library `MUMPS_freefem.dylib` (MacOSX), `MUMPS_freefem.so` (Unix) or `MUMPS_freefem.dll` (Windows).

This is done in typing `:::freefem load "MUMPS_FreeFem"` in the `.edp` file. We give now the two methods to give the option of MUMPS solver in __`FreeFem++`__.

* __Solver parameters is defined in .edp file:__
	In this method, we need to give the parameters `:::freefem lparams` and `:::freefem dparams`. These parameters are defined for MUMPS by :

	`:::freefem lparams[0] = SYM`,<br>
	`:::freefem lparams[1] = PAR`,<br>
	$\forall i$ = 1,...,40, `:::freefem lparams[i+1] = ICNTL(i)`<br>
	$\forall i$ = 1,...,15, `:::freefem dparams[i-1] = CNTL(i)`

* __Reading solver parameters on a file:__

	The structure of data file for MUMPS in __`FreeFem++`__ is : first line parameter `SYM` and second line parameter `PAR` and in the following line the different value of vectors `ICNTL` and `CNTL`. An example of this parameter file is given in `:::freefem ffmumpsfileparam.txt`.

	```freefem
	0 /* SYM :: 0 for non symmetric matrix, 1 for symmetric definite positive matrix and 2 general symmetric matrix*/
	1 /* PAR :: 0 host not working during factorization and solves steps, 1 host working during factorization and solves steps*/
	-1 /* ICNTL(1) :: output stream for error message */
	-1 /* ICNTL(2) :: output for diagnostic printing, statics and warning message */
	-1 /* ICNTL(3) :: for global information */
	0 /* ICNTL(4) :: Level of printing for error, warning and diagnostic message */
	0 /* ICNTL(5) :: matrix format : 0 assembled format, 1 elemental format. */
	7 /* ICNTL(6) :: control option for permuting and/or scaling the matrix in analysis phase */
	3 /* ICNTL(7) :: pivot order strategy : AMD, AMF, metis, pord scotch*/
	77 /* ICNTL(8) :: Row and Column scaling strategy */
	1 /* ICNTL(9) :: 0 solve Ax = b, 1 solve the transposed system A^t x = b : parameter is not considered in the current release of freefem++*/
	0 /* ICNTL(10) :: number of steps of iterative refinement */
	0 /* ICNTL(11) :: statics related to linear system depending on ICNTL(9) */
	1 /* ICNTL(12) :: constrained ordering strategy for general symmetric matrix */
	0 /* ICNTL(13) :: method to control splitting of the root frontal matrix */
	20 /* ICNTL(14) :: percentage increase in the estimated working space (default 20\%)*/
	0 /* ICNTL(15) :: not used in this release of MUMPS */
	0 /* ICNTL(16) :: not used in this release of MUMPS */
	0 /* ICNTL(17) :: not used in this release of MUMPS */
	3 /* ICNTL(18) :: method for given : matrix pattern and matrix entries : */
	0 /* ICNTL(19) :: method to return the Schur complement matrix */
	0 /* ICNTL(20) :: right hand side form ( 0 dense form, 1 sparse form) : parameter will be set to 0 for freefem++ */
	0 /* ICNTL(21) :: 0, 1 kept distributed solution : parameter is not considered in the current release of freefem++ */
	0 /* ICNTL(22) :: controls the in-core/out-of-core (OOC) facility */
	0 /* ICNTL(23) :: maximum size of the working memory in Megabyte than MUMPS can allocate per working processor */
	0 /* ICNTL(24) :: control the detection of null pivot */
	0 /* ICNTL(25) :: control the computation of a null space basis */
	0 /* ICNTL(26) :: This parameter is only significant with Schur option (ICNTL(19) not zero). : parameter is not considered in the current release of freefem++ */
	-8 /* ICNTL(27) (Experimental parameter subject to change in next release of MUMPS) :: control the blocking factor for multiple righthand side during the solution phase : parameter is not considered in the current release of freefem++ */
	0 /* ICNTL(28) :: not used in this release of MUMPS*/
	0 /* ICNTL(29) :: not used in this release of MUMPS*/
	0 /* ICNTL(30) :: not used in this release of MUMPS*/
	0 /* ICNTL(31) :: not used in this release of MUMPS*/
	0 /* ICNTL(32) :: not used in this release of MUMPS*/
	0 /* ICNTL(33) :: not used in this release of MUMPS*/
	0 /* ICNTL(34) :: not used in this release of MUMPS*/
	0 /* ICNTL(35) :: not used in this release of MUMPS*/
	0 /* ICNTL(36) :: not used in this release of MUMPS*/
	0 /* ICNTL(37) :: not used in this release of MUMPS*/
	0 /* ICNTL(38) :: not used in this release of MUMPS*/
	1 /* ICNTL(39) :: not used in this release of MUMPS*/
	0 /* ICNTL(40) :: not used in this release of MUMPS*/
	0.01 /* CNTL(1) :: relative threshold for numerical pivoting */
	1e-8 /* CNTL(2) :: stopping criteria for iterative refinement */
	-1 /* CNTL(3) :: threshold for null pivot detection */
	-1 /* CNTL(4) :: determine the threshold for partial pivoting */
	0.0 /* CNTL(5) :: fixation for null pivots */
	0 /* CNTL(6) :: not used in this release of MUMPS */
	0 /* CNTL(7) :: not used in this release of MUMPS */
	0 /* CNTL(8) :: not used in this release of MUMPS */
	0 /* CNTL(9) :: not used in this release of MUMPS */
	0 /* CNTL(10) :: not used in this release of MUMPS */
	0 /* CNTL(11) :: not used in this release of MUMPS */
	0 /* CNTL(12) :: not used in this release of MUMPS */
	0 /* CNTL(13) :: not used in this release of MUMPS */
	0 /* CNTL(14) :: not used in this release of MUMPS */
	0 /* CNTL(15) :: not used in this release of MUMPS */
	```

If no solver parameter is given, we used default option of MUMPS solver.

!!!question "MUMPS example"
	A simple example of calling MUMPS in __`FreeFem++`__ with this two methods is given in the [Test solver MUMPS example](../examples/#solver-mumps).

### SuperLU distributed solver

The package [SuperLU_DIST](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) solves linear systems using LU factorization. It is a free scientific library under BSD license.

This library provides functions to handle square or rectangular matrix in real and complex arithmetics. The method implemented in SuperLU_DIST is a supernodal method. New release of this package includes a parallel symbolic factorization. This scientific library is written in C and MPI for communications.

__SuperLU_DIST parameters:__

We describe now some parameters of SuperLU_DIST. The SuperLU_DIST library use a 2D-logical process group. This process grid is specified by $nprow$ (process row) and $npcol$ (process column) such that $N_{p} = nprow \: npcol$ where $N_{p}$ is the number of all process allocated for SuperLU_DIST.

The input matrix parameters is controlled by "matrix= " in `:::freefem sparams` for internal parameter or in the third line of parameters file. The different value are

* `:::freefem matrix=assembled` global matrix are available on all process
* `:::freefem matrix=distributedglobal` The global matrix is distributed among all the process
* `:::freefem matrix=distributed` The input matrix is distributed (not yet implemented)

The option arguments of SuperLU_DIST are described in the section Users-callable routine of the [SuperLU users' guide](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ug.pdf).

The parameter `Fact` and `TRANS` are specified in __`FreeFem++`__ interfaces to SuperLU_DIST during the different steps. For this reason, the value given by the user for this option is not considered.

The factorization LU is calculated in SuperLU_DIST on the matrix $A_p$.
$$
A_{p} = P_{c} \: P_r \: D_r \: A \: D_{c} \: P_{c}^{t}
$$
where $P_c$ and $P_r$ is the row and column permutation matrix respectively, $D_r$ and $D_c$ are diagonal matrix for respectively row and column scaling.

The option argument `RowPerm` (resp. `ColPerm`) control the row (resp. column) permutation matrix. $D_r$ and $D_c$ is controlled by the parameter `DiagScale`.

The parameter `:::freefem permr`, `:::freefem permc`, `:::freefem scaler`, and `:::freefem scalec` in __`FreeFem++`__ is provided to give row permutation, column permutation, row scaling and column scaling of the user respectively.

The other parameters for LU factorization are `ParSymFact` and `ReplaceTinyPivot`. The parallel symbolic factorization works only on a power of two processes and need the `ParMetis` ordering. The default option argument of SuperLU_DIST are given in the file `ffsuperlu_dist_fileparam.txt`.

__Calling SuperLU_DIST in __`FreeFem++`____

To call SuperLU_DIST in __`FreeFem++`__, we need to load the library dynamic correspond to interface. This done by the following line `:::freefem load "real_superlu _DIST_FreeFem"` (resp. `:::freefem load "complex_superlu_DIST_FreeFem"`) for real (resp. complex) arithmetics in the file `.edp`.

__Solver parameters is defined in `.edp` file:__

To call SuperLU_DIST with internal parameter, we used the parameters `sparams`. The value of parameters of SuperLU_DIST in `sparams` are defined by :

* `nprow=1`,
* `npcol=1`,
* `matrix= distributedgloba`,
* `Fact= DOFACT`,
* `Equil=NO`,
* `ParSymbFact=NO`,
* `ColPerm= MMD_AT_PLUS_A`,
* `RowPerm= LargeDiag`,
* `DiagPivotThresh=1.0`,
* `IterRefine=DOUBLE`,
* `Trans=NOTRANS`,
* `ReplaceTinyPivot=NO`,
* `SolveInitialized=NO`,
* `PrintStat=NO`,
* `DiagScale=NOEQUIL`

This value correspond to the parameter in the file `ffsuperlu_dist_fileparam.txt`. If one parameter is not specified by the user, we take the default value of SuperLU_DIST.

__Reading solver parameters on a file:__
The structure of data file for SuperLU_DIST in __`FreeFem++`__ is given in the file `ffsuperlu_dist_fileparam.txt` (default value of the __`FreeFem++`__ interface).

```freefem
1 /* nprow : integer value */
1 /* npcol : integer value */
distributedglobal /* matrix input : assembled, distributedglobal, distributed */
DOFACT /* Fact : DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED */
NO /* Equil : NO, YES */
NO /* ParSymbFact : NO, YES */
MMD_AT_PLUS_A /* ColPerm : NATURAL, MMD_AT_PLUS_A, MMD_ATA, METIS_AT_PLUS_A, PARMETIS, MY_PERMC */
LargeDiag /* RowPerm : NOROWPERM, LargeDiag, MY_PERMR */
1.0 /* DiagPivotThresh : real value */
DOUBLE /* IterRefine : NOREFINE, SINGLE, DOUBLE, EXTRA */
NOTRANS /* Trans : NOTRANS, TRANS, CONJ*/
NO /* ReplaceTinyPivot : NO, YES*/
NO /* SolveInitialized : NO, YES*/
NO /* RefineInitialized : NO, YES*/
NO /* PrintStat : NO, YES*/
NOEQUIL /* DiagScale : NOEQUIL, ROW, COL, BOTH*/
```

If no solver parameter is given, we used default option of SuperLU_DIST solver.

!!!question "Example"
	A simple example of calling SuperLU_DIST in __`FreeFem++`__ with this two methods is given in the [Solver superLU_DIST example](../examples/#solver-superlu_dist).

### PaStiX solver

[PaStiX](http://pastix.gforge.inria.fr/files/README-txt.html) (Parallel Sparse matrix package) is a free scientific library under CECILL-C license. This package solves sparse linear system with a direct and block ILU(k) iterative methods. This solver can be applied to a real or complex matrix with a symmetric pattern.

__PaStiX parameters:__

The input `:::freefem matrix` parameter of __`FreeFem++`__ depend on PaStiX interface. `:::freefem matrix = assembled` for non distributed matrix. It is the same parameter for SuperLU_DIST.

There are four parameters in PaStiX : `iparm`, `dparm`, `perm` and `invp`. These parameters are respectively the integer parameters (vector of size 64), real parameters (vector of size 64), permutation matrix and inverse permutation matrix respectively. `iparm` and `dparm` vectors are described in [PaStiX RefCard](https://gforge.inria.fr/docman/?group_id=186&view=listfile&dirid=246).

The parameters `:::freefem permr` and `:::freefem permc` in __`FreeFem++`__ are provided to give permutation matrix and inverse permutation matrix of the user respectively.

__Solver parameters defined in `.edp` file:__

To call PaStiX in __`FreeFem++`__ in this case, we need to specify the parameters `:::freefem lparams` and `:::freefem dparams`. These parameters are defined by :

$\forall i$ = 0,... ,63, `:::freefem lparams[i] = iparm[i]`.

$\forall i$ = 0,... ,63, `:::freefem dparams[i] = dparm[i]`.


__Reading solver parameters on a file:__

The structure of data file for PaStiX parameters in __`FreeFem++`__ is : first line structure parameters of the matrix and in the following line the value of vectors `iparm` and `dparm` in this order.

```freefem
assembled /* matrix input :: assembled, distributed global and distributed */
iparm[0]
iparm[1]
...
...
iparm[63]
dparm[0]
dparm[1]
...
...
dparm[63]
```

An example of this file parameter is given in `ffpastix_iparm_dparm.txt` with a description of these parameters. This file is obtained with the example file `iparm.txt` and `dparm.txt` including in the PaStiX package.

If no solver parameter is given, we use the default option of PaStiX solver.

!!!question "Example"
	A simple example of calling PaStiX in __`FreeFem++`__ with this two methods is given in the [Solver PaStiX example](../examples/#solver-pastix).

In [Table 3](#Tab3), we recall the different matrix considering in the different direct solvers.

<table>
	<thead>
		<tr>
			<th colspan="7"><a name="Tab3">Table 3</a>: Type of matrix used by the different direct sparse solver</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td rowspan="2" style="vertical-align: bottom !important; font-weight: bold">direct solver</td>
			<td colspan="3" align="center" style="font-weight: bold">square matrix</td>
			<td colspan="3" align="center" style="font-weight: bold">rectangular matrix</td>
		</tr>
		<tr>
			<td style="font-weight: bold">sym</td>
			<td align="center" style="font-weight: bold">sym pattern</td>
			<td align="center" style="font-weight: bold">unsym</td>
			<td align="center" style="font-weight: bold">sym</td>
			<td align="center" style="font-weight: bold">sym pattern</td>
			<td align="center" style="font-weight: bold">unsym</td>
		</tr>
		<tr>
			<td>SuperLU_DIST</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
		</tr>
		<tr>
			<td>MUMPS</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">no</td>
			<td align="center">no</td>
			<td align="center">no</td>
		</tr>
		<tr>
			<td>pastix</td>
			<td align="center">yes</td>
			<td align="center">yes</td>
			<td align="center">no</td>
			<td align="center">no</td>
			<td align="center">no</td>
			<td align="center">no</td>
		</tr>
	</tbody>
</table>


## Parallel sparse iterative solver

Concerning iterative solvers, we have chosen [pARMS](http://www-users.cs.umn.edu/~saad/software/pARMS/), [HIPS](http://hips.gforge.inria.fr/) and [Hypre](https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods).

Each software implements a different type of parallel preconditioner.

So, pARMS implements algebraic domain decomposition preconditioner type such as additive Schwartz [CAI1989](#CAI1989) and interface method; while HIPS implement hierarchical incomplete factorization and finally HYPRE implements multilevel preconditioner are AMG(Algebraic MultiGrid) and parallel approximated inverse.

To use one of these programs in __`FreeFem++`__, you have to install it independently of __`FreeFem++`__. It is also necessary to install the MPI communication library which is essential for communication between the processors and, in some cases, software partitioning graphs like [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) or [Scotch](http://www.labri.fr/perso/pelegrin/scotch/).

All this preconditioners are used with Krylov subspace methods accelerators.

Krylov subspace methods are iterative methods which consist in finding a solution $x$ of linear system $Ax=b$ inside the affine space $x_0+K_m$ by imposing that $b-Ax \bot \mathcal{L}_m$, where $K_m$ is Krylov subspace of dimension $m$ defined by $K_m=\{r_0, Ar_0, A^2r_0,...,A^{m-1}r_0\}$ and $\mathcal{L}_m$ is another subspace of dimension $m$ which depends on type of Krylov subspace. For example in GMRES, $\mathcal{L}_m=AK_m$.

We realized an interface which is easy to use, so that the call of these different softwares in __`FreeFem++`__ is done in the same way. You just have to load the solver and then specify the parameters to apply to the specific solvers. In the rest of this chapter, when we talk about Krylov subspace methods we mean one among GMRES, CG and BICGSTAB.

### pARMS solver

[pARMS](http://www-users.cs.umn.edu/~saad/software/pARMS/) (parallel Algebraic Multilevel Solver) is a software developed by Youssef Saad and al at University of Minnesota.

This software is specialized in the resolution of large sparse non symmetric linear systems of equation. Solvers developed in pARMS are of type "Krylov's subspace".

It consists of variants of GMRES like FGMRES (Flexible GMRES), DGMRES (Deflated GMRES) [SAAD2003](#SAAD2003) and BICGSTAB. pARMS also implements parallel preconditioner like RAS (Restricted Additive Schwarz) [CAI1989](#CAI1989) and Schur Complement type preconditioner.

All these parallel preconditioners are based on the principle of domain decomposition. Thus, the matrix $A$ is partitioned into sub matrices $A_i$($i=1,...,p$) where p represents the number of partitions one needs. The union of $A_i$ forms the original matrix. The solution of the overall system is obtained by solving the local systems on $A_i$ (see [SMITH1996](#SMITH1996)). Therefore, a distinction is made between iterations on $A$ and the local iterations on $A_i$.

To solve the local problem on $A_i$ there are several preconditioners as __ilut__ (Incomplete LU with threshold), __iluk__ (Incomplete LU with level of fill in) and __ARMS__ (Algebraic Recursive Multilevel Solver).

!!!question "Default parameters"
	```freefem
	load "parms_FreeFem" //Tell FreeFem that you will use pARMS

	// Mesh
	border C(t=0, 2*pi){x=cos(t); y=sin(t); label=1;}
	mesh Th = buildmesh (C(50));

	// Fespace
	fespace Vh(Th, P2);
	Vh u, v;

	// Function
	func f= x*y;

	// Problem
	problem Poisson (u, v, solver=sparsesolver)
		= int2d(Th)(
			  dx(u)*dx(v)
			+ dy(u)*dy(v)
		)
		+ int2d(Th)(
			- f*v
		)
		+ on(1, u=0)
		;

	// Solve
	real cpu = clock();
	Poisson;
	cout << " CPU time = " << clock()-cpu << endl;

	// Plot
	plot(u);
	```

	In line 1, the pARMS dynamic library is loaded with interface __`FreeFem++`__. After this, in line 15 we specify that the bilinear form will be solved by the last sparse linear solver load in memory which, in this case, is pARMS.

	The parameters used in pARMS in this case are the default one since the user does not have to provide any parameter.

Here are some default parameters:

* `solver=FGMRES`,
* `Krylov dimension=30`,
* `Maximum of Krylov=1000`,
* `Tolerance for convergence=$1e-08$`(see book [SAAD2003](#SAAD2003) to understand all this parameters),
* `preconditionner=Restricted Additif Schwarz` [CAI1989](#CAI1989),
* `Inner Krylov dimension=5`,
* `Maximum of inner Krylov dimension=5`,
* `Inner preconditionner=ILUK`.

To specify the parameters to apply to the solver, the user can either give an integer vector for __integer parameters__ and real vectors for __real parameters__ or provide a __file__ which contains those parameters.

!!!question "User specifies parameters inside two vectors"
	Lets us consider Navier-Stokes example. In this example we solve linear systems coming from discretization of Navier-Stokes equations with pARMS. Parameters of solver is specified by user.

	```freefem
	load "parms_FreeFem"

	// Parameters
	real nu = 1.;
	int[int] iparm(16);
	real[int] dparm(6);
	for (int ii = 0; ii < 16; ii++)
		iparm[ii] = -1;
	for (int ii = 0; ii < 6; ii++)
		dparm[ii] = -1.0;
	iparm[0]=0;

	// Mesh
	mesh Th = square(10, 10);
	int[int] wall = [1, 3];
	int inlet = 4;

	// Fespace
	fespace Vh(Th, [P2, P2, P1]);

	// Function
	func uc = 1.;

	varf Stokes ([u, v, p], [ush, vsh, psh], solver=sparsesolver)
		= int2d(Th)(
			  nu*(
				  dx(u)*dx(ush)
				+ dy(u)*dy(ush)
				+ dx(v)*dx(vsh)
				+ dy(v)*dy(vsh)
			)
			- p*psh*(1.e-6)
			- p*(dx(ush) + dy(vsh))
			- (dx(u) + dy(v))*psh
		)
		+ on(wall, wall, u=0., v=0.)
		+ on(inlet, u=uc, v=0)
		;

	matrix AA = Stokes(Vh, Vh);
	set(AA, solver=sparsesolver, lparams=iparm, dparams=dparm); //set pARMS as linear solver
	real[int] bb = Stokes(0, Vh);
	real[int] sol(AA.n);
	sol = AA^-1 * bb;
	```

	We need two vectors to specify the parameters of the linear solver. In line 5-6 of the example, we have declared these vectors(`:::freefem int[int] iparm(16); real[int] dparm(6);`). In line 7-10 we have initialized these vectors by negative values.

	We do this because all parameters values in pARMS are positive and if you do not change the negative values of one entry of this vector, the default value will be set.

	In [table 4](#Tab4) and [table 5](#Tab5), we have the meaning of different entries of these vectors.

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab4">Table 4</a>: Meaning of `:::freefem lparams` corresponding variables</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td style="font-weight: bold">Entries of `:::freefem iparm`</td>
				<td style="font-weight: bold">Significations of each entries</td>
			</tr>
			<tr>
				<td>`iparm[0]`</td>
				<td>Krylov subspace methods.<br>Different values for this parameters are specify on [table 7](#Tab6)</td>
			</tr>
			<tr>
				<td>`iparm[1]`</td>
				<td>Preconditionner.<br>Different preconditionners for this parameters are specify on [table 7](#Tab7)</td>
			</tr>
			<tr>
				<td>`iparm[2]`</td>
				<td>Krylov subspace dimension in outer iteration: default value 30</td>
			</tr>
			<tr>
				<td>`iparm[3]`</td>
				<td>Maximum of iterations in outer iteration: default value 1000</td>
			</tr>
			<tr>
				<td>`iparm[4]`</td>
				<td>Number of level in arms when used.</td>
			</tr>
			<tr>
				<td>`iparm[5]`</td>
				<td>Krylov subspace dimension in inner iteration: default value 3</td>
			</tr>
			<tr>
				<td>`iparm[6]`</td>
				<td>Maximum of iterations in inner iteration: default value 3</td>
			</tr>
			<tr>
				<td>`iparm[7]`</td>
				<td>Symmetric(=1 for symmetric) or unsymmetric matrix:<br>default value 0(unsymmetric matrix)</td>
			</tr>
			<tr>
				<td>`iparm[8]`</td>
				<td>Overlap size between different subdomain: default value 0(no overlap)</td>
			</tr>
			<tr>
				<td>`iparm[9]`</td>
				<td>Scale the input matrix or not: Default value 1 (Matrix should bescale)</td>
			</tr>
			<tr>
				<td>`iparm[10]`</td>
				<td>Block size in arms when used: default value 20</td>
			</tr>
			<tr>
				<td>`iparm[11]`</td>
				<td>lfil0 (ilut, iluk, and arms) : default value 20</td>
			</tr>
			<tr>
				<td>`iparm[12]`</td>
				<td>lfil for Schur complement const : default value 20</td>
			</tr>
			<tr>
				<td>`iparm[13]`</td>
				<td>lfil for Schur complement const : default value 20</td>
			</tr>
			<tr>
				<td>`iparm[14]`</td>
				<td>Multicoloring or not in ILU when used : default value 1</td>
			</tr>
			<tr>
				<td>`iparm[15]`</td>
				<td>Inner iteration : default value 0</td>
			</tr>
			<tr>
				<td>`iparm[16]`</td>
				<td>Print message when solving:default 0 (no message print).<br>0: no message is print,<br>1: Convergence informations like number of iteration and residual ,<br>2: Timing for a different step like preconditioner<br>3 : Print all informations.</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab5">Table 5</a>: Significations of `:::freefem dparams` corresponding variables</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td style="font-weight: bold">Entries of `:::freefem dparm`</td>
				<td style="font-weight: bold">Significations of each entries</td>
			</tr>
			<tr>
				<td>`dparm[0]`</td>
				<td>precision for outer iteration : default value 1e-08</td>
			</tr>
			<tr>
				<td>`dparm[1]`</td>
				<td>precision for inner iteration: default value 1e-2</td>
			</tr>
			<tr>
				<td>`dparm[2]`</td>
				<td>tolerance used for diagonal domain: : default value 0.1</td>
			</tr>
			<tr>
				<td>`dparm[3]`</td>
				<td>drop tolerance droptol0 (ilut, iluk, and arms) : default value 1e-2</td>
			</tr>
			<tr>
				<td>`dparm[4]`</td>
				<td>droptol for Schur complement const: default value 1e-2</td>
			</tr>
			<tr>
				<td>`dparm[5]`</td>
				<td>droptol for Schur complement const: default value 1e-2</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab6">Table 6</a>: Krylov Solvers in pARMS</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td style="font-weight: bold">Values of `iparm[0]`</td>
				<td style="font-weight: bold">Krylov subspace methods</td>
			</tr>
			<tr>
				<td>0</td>
				<td>FGMRES (Flexible GMRES)</td>
			</tr>
			<tr>
				<td>1</td>
				<td>DGMRES (Deflated GMRES)</td>
			</tr>
			<tr>
				<td>2</td>
				<td>BICGSTAB</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Tab7">Table 7</a>: Preconditionners in pARMS</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td style="font-weight: bold">Values of `iparm[1]`</td>
				<td style="font-weight: bold">Preconditionners type</td>
			</tr>
			<tr>
				<td>0</td>
				<td>additive Schwartz preconditioner with ilu0 as local preconditioner</td>
			</tr>
			<tr>
				<td>1</td>
				<td>additive Schwartz preconditioner with iluk as local preconditioner</td>
			</tr>
			<tr>
				<td>2</td>
				<td>additive Schwartz preconditioner with ilut as local preconditioner</td>
			</tr>
			<tr>
				<td>3</td>
				<td>additive Schwartz preconditioner with arms as local preconditioner</td>
			</tr>
			<tr>
				<td>4</td>
				<td>Left Schur complement preconditioner with ilu0 as local preconditioner</td>
			</tr>
			<tr>
				<td>5</td>
				<td>Left Schur complement preconditioner with ilut as local preconditioner</td>
			</tr>
			<tr>
				<td>6</td>
				<td>Left Schur complement preconditioner with iluk as local preconditioner</td>
			</tr>
			<tr>
				<td>7</td>
				<td>Left Schur complement preconditioner with arms as local preconditioner</td>
			</tr>
			<tr>
				<td>8</td>
				<td>Right Schur complement preconditioner with ilu0 as local preconditioner</td>
			</tr>
			<tr>
				<td>9</td>
				<td>Right Schur complement preconditioner with ilut as local preconditioner</td>
			</tr>
			<tr>
				<td>10</td>
				<td>Right Schur complement preconditioner with iluk as local preconditioner</td>
			</tr>
			<tr>
				<td>11</td>
				<td>Right Schur complement preconditioner with arms as local preconditioner</td>
			</tr>
			<tr>
				<td>12</td>
				<td>sch_gilu0, Schur complement preconditioner with global ilu0</td>
			</tr>
			<tr>
				<td>13</td>
				<td>SchurSymmetric GS preconditioner</td>
			</tr>
		</tbody>
	</table>

	We run this example on a cluster paradent of Grid5000 and report results in [table 8](#Tab8).

	<table>
		<thead>
			<tr>
				<th colspan="5"><a name="Tab8">Table 8</a>: Convergence and time for solving linear system</th>
			</tr>
		</thead>
		<tbody>
			<tr>
				<td colspan="2" align="center" style="font-weight: bold">n=471281</td>
				<td colspan="2" align="center" style="font-weight: bold">nnz=$13\times10^6$</td>
				<td colspan="2" align="center" style="font-weight: bold">Te=571,29</td>
			</tr>
			<tr>
				<td rowspan="2" align="center" style="vertical-align: bottom !important">np</td>
				<td colspan="2" align="center">add(iluk)</td>
				<td colspan="2" align="center">schur(iluk)</td>
			</tr>
			<tr>
				<td align="center">`nit`</td>
				<td align="center">`time`</td>
				<td align="center">`nit`</td>
				<td align="center">`time`</td>
			</tr>
			<tr>
				<td align="center">4</td>
				<td align="center">230</td>
				<td align="center">637.57</td>
				<td align="center">21</td>
				<td align="center">557.8</td>
			</tr>
			<tr>
				<td align="center">8</td>
				<td align="center">240</td>
				<td align="center">364.12</td>
				<td align="center">22</td>
				<td align="center">302.25</td>
			</tr>
			<tr>
				<td align="center">16</td>
				<td align="center">247</td>
				<td align="center">212.07</td>
				<td align="center">24</td>
				<td align="center">167.5</td>
			</tr>
			<tr>
				<td align="center">32</td>
				<td align="center">261</td>
				<td align="center">111.16</td>
				<td align="center">25</td>
				<td align="center">81.5</td>
			</tr>
		</tbody>
	</table>

	<table>
		<thead>
			<tr>
				<th colspan="2"><a name="Ta9">Table 9</a>: Legend of [table 8](#Tab8)</th>
			</tr>
		</thead>
		<tbody>
		<tr>
			<td>`n`</td>
			<td>matrix size</td>
		</tr>
		<tr>
			<td>`nnz`</td>
			<td>number of non null entries inside matrix</td>
		</tr>
		<tr>
			<td>`nit`</td>
			<td>number of iteration for convergence</td>
		</tr>
		<tr>
			<td>`time`</td>
			<td>Time for convergence</td>
		</tr>
		<tr>
			<td>`Te`</td>
			<td>Time for constructing finite element matrix</td>
		</tr>
		<tr>
			<td>`np`</td>
			<td>number of processor</td>
		</tr>
	</table>

	In this example, we fix the matrix size (in term of finite element, we fix the mesh) and increase the number of processors used to solve the linear system. We saw that, when the number of processors increases, the time for solving the linear equation decreases, even if the number of iteration increases. This proves that, using pARMS as solver of linear systems coming from discretization of partial differential equation in __`FreeFem++`__ can decrease drastically the total time of simulation.

## References

<a name="CAI1989">[CAI1989]</a> CAI, Xiao-Chuan. Some domain decomposition algorithms for nonselfadjoint elliptic and parabolic partial differential equations. 1989.

<a name="SAAD2003">[SAAD2003]</a> SAAD, Yousef. Iterative methods for sparse linear systems. siam, 2003.

<a name="SMITH1996">[SMITH1996]</a> SMITH, B. P. Bj rstad and W. Gropp, Domain Decomposition. 1996.
