# Parallel sparse solvers

Parallel sparse solvers use several processors to solve linear systems of equation. Like sequential,
parallel linear solvers can be direct or iterative. In FreeFem++ both are available.

## Using parallel sparse solvers in FreeFem++

We recall that the `:::freefem solver` parameters are defined in the following commands: `:::freefem solve`, `:::freefem problem`, `:::freefem set` (setting parameter of a matrix) and in the construction of the matrix corresponding to a bilinear form. In these commands, the parameter `:::freefem solver` must be set to `:::freefem sparsesolver` for parallel sparse solver. We have added specify parameters to these command lines for parallel sparse solvers. These are

* `:::freefem lparams` : vector of integer parameters (l is for the c++ type long)
* `:::freefem dparams` : vector of real parameters
* `:::freefem sparams` : string parameters
* `:::freefem datafilename` : name of the file which contains solver parameters

The following four parameters are only for direct solvers and are vectors. These parameters allow the user to preprocess the matrix (see the section on sparse direct solver above for more information).

* `:::freefem permr` : row permutation (integer vector)
* `:::freefem permc` : column permutation or inverse row permutation (integer vector)
* `:::freefem scaler` : row scaling (real vector)
* `:::freefem scalec` : column scaling (real vector)


There are two possibilities to control solver parameters. The first method defines parameters with `:::freefem lparams`, `:::freefem dparams` and `:::freefem sparams` in .edp file. The second one reads the solver parameters from a data file. The name of this file is specified by `:::freefem datafilename`. If `:::freefem lparams`, `:::freefem dparams`, `:::freefem sparams` or `:::freefem datafilename` is not provided by the user, the solver's default value is used.

To use parallel solver in FreeFem++, we need to load the dynamic library corresponding to this solver.
For example to use MUMPS solver as parallel solver in FreeFem, write in the .edp file `:::freefem load "MUMPS\_FreeFem"`.

If the libraries are not loaded, the default sparse solver will be loaded (default sparse solver is UMFPACK). The table 11.1 \ref{lib.sparse.solver} $\codered$ gives this new value for the different libraries.

<table>
	<thead>
		<tr>
			<th colspan="3">Table 11.1 : Default sparse solver for real and complex arithmetics when we load a parallel sparse solver library</th>
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

We also add functions (see Table 11.2 \ref{func.sparse.solver} $\codered$) with no parameter to change the default sparse solver in the .edp file. To use these functions, we need to load the library corresponding to the solver. An example of using different parallel sparse solvers for the same problem is given in testdirectsolvers.edp (directory example$++-$mpi $\codered$).

<table>
	<thead>
		<tr>
			<th colspan="3">Table 11.2 : Functions that allow to change the default sparse solver for real and complex arithmetics and the result of these functions</th>
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

__Example testdirectsolvers.edp__

```freefem
load "../src/solver/MUMPS_FreeFem"
// default solver : real-> MUMPS, complex -> MUMPS
load "../src/solver/real_SuperLU_DIST_FreeFem"
// default solver : real-> SuperLU_DIST, complex -> MUMPS
load "../src/solver/real_pastix_FreeFem"
// default solver : real-> pastix, complex -> MUMPS

// solving with pastix
{
 matrix A =
 [[ 1, 2, 2, 1, 1],
 [ 2, 12, 0, 10 , 10],
 [ 2, 0, 1, 0, 2],
 [ 1, 10, 0, 22, 0.],
 [ 1, 10, 2, 0., 22]];

 real[int] xx = [ 1,32,45,7,2], x(5), b(5), di(5);
 b=A*xx;
 cout << "b=" << b << endl;
 cout << "xx=" << xx << endl;

 set(A,solver=sparsesolver,datafilename="ffpastix_iparm_dparm.txt");
 cout << "solving solution" << endl;
 x = A^-1*b;
 cout << "b=" << b << endl;
 cout << "x=" << endl; cout << x << endl;
 di = xx-x;
 if(mpirank==0){
 cout << "x-xx="<< endl; cout << "Linf "<< di.linfty << " L2 " << di.l2 << endl;
 }
}

// solving with SuperLU_DIST
realdefaulttoSuperLUdist();
// default solver : real-> SuperLU_DIST, complex -> MUMPS
{
 matrix A =
 [[ 1, 2, 2, 1, 1],
 [ 2, 12, 0, 10 , 10],
 [ 2, 0, 1, 0, 2],
 [ 1, 10, 0, 22, 0.],
 [ 1, 10, 2, 0., 22]];

 real[int] xx = [ 1,32,45,7,2], x(5), b(5), di(5);
 b=A*xx;
 cout << "b=" << b << endl;
 cout << "xx=" << xx << endl;

 set(A,solver=sparsesolver,datafilename="ffsuperlu_dist_fileparam.txt");
 cout << "solving solution" << endl;
 x = A^-1*b;
 cout << "b=" << b << endl;
 cout << "x=" << endl; cout << x << endl;
 di = xx-x;
 if(mpirank==0){
 cout << "x-xx="<< endl; cout << "Linf "<< di.linfty << " L2 " << di.l2 << endl;
 }
}

// solving with MUMPS
defaulttoMUMPS();
// default solver : real-> MUMPS, complex -> MUMPS
{
 matrix A =
 [[ 1, 2, 2, 1, 1],
 [ 2, 12, 0, 10 , 10],
 [ 2, 0, 1, 0, 2],
 [ 1, 10, 0, 22, 0.],
 [ 1, 10, 2, 0., 22]];

 real[int] xx = [ 1,32,45,7,2], x(5), b(5), di(5);
 b=A*xx;
 cout << "b=" << b << endl;
 cout << "xx=" << xx << endl;

 set(A,solver=sparsesolver,datafilename="ffmumps_fileparam.txt");
 cout << "solving solution" << endl;
 x = A^-1*b;
 cout << "b=" << b << endl;
 cout << "x=" << endl; cout << x << endl;
 di = xx-x;
 if(mpirank==0){
 cout << "x-xx="<< endl; cout << "Linf "<< di.linfty << " L2 " << di.l2 << endl;
 }
}
```

$\codered$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A faire :: \\
% M\"{\i}?`ï¿½thode Freefem++ r\"{\i}?`ï¿½solution de solveur : comment les options des parallels solveurs dans freefem++
%d\"{\i}?`ï¿½finir: lparams, dparams, ..... lecture d'un fichier.
%comment utiliser tel ou tel solveur dans Freefem++
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Sparse direct solver

In this section, we present the sparse direct solvers interfaced with FreeFem++.

### MUMPS solver

MUltifrontal Massively Parallel Solver (MUMPS) is a free library \cite{mumpspubl1,mumpspubl2,mumpspubl3} $\codered$. This package solves linear system of the form $A \: x = b$ where $A$ is a square sparse matrix with a direct method. The square matrix considered in MUMPS can be either unsymmetric, symmetric positive definite or general symmetric.
The method implemented in MUMPS is a direct method based on a multifrontal approach \cite{mumpspubl1} $\codered$. It constructs a direct factorization $A \:= \: L\:U$, $A\: = \: L^t \: D \: L$ depending of the symmetry of the matrix $A$. MUMPS uses the following libraries : BLAS\cite{blas1,blas2} $\codered$, BLACS and ScaLAPACK\cite{scalapackuserguide} $\codered$.

!!! note
	MUMPS does not solve linear system with a rectangular matrix.



__Installation of MUMPS__

To used MUMPS in FreeFem++, you have to install the MUMPS package into your computer.
MUMPS is written in Fortran 90. The parallel version is constructed using MPI \cite{mpi} $\codered$ for message passing and BLAS \cite{blas1,blas2} $\codered$, BLACS and ScaLAPACK\cite{scalapackuserguide} $\codered$. Therefore, a fortran compiler is needed, and MPI, BLAS, BLACS and ScaLAPACK. An installation procedure to obtain this package is given in the file README\_COMPILE in the directory src/solver $\codered$ of FreeFem++.

__Creating Library of MUMPS interface for FreeFem++:__

The MUMPS interface for FreeFem++ is given in file MUMPS\_freefem.cpp (directory src/solver/ $\codered$).
This interface works with the release 3.8.3 and 3.8.4 of MUMPS. To used MUMPS in FreeFem++, we need the library corresponding to this interface.
A description to obtain this library is given in the file README\_COMPILE in the directory src/solver $\codered$ of FreeFem++. We recall here the procedure. Go to the directory src/solver in FreeFem++ package. Edit the file makefile-sparsesolver.inc to yours system: comment Section 1, comment line corresponding to libraries BLAS, BLACS, ScaLAPACK, Metis, scotch in Section 2 and comment in Section 3 $\codered$ the paragraph corresponding to MUMPS solver. And then type `:::bash make mumps` in a terminal window.

Now we give a short description of MUMPS parameters before describing the method to call MUMPS in FreeFem++.

__MUMPS parameters:__

There are four input parameters in MUMPS (see \cite{mumpsuserguide} $\codered$). Two integers SYM and PAR, a vector of integer of size 40 INCTL and a vector of real of size 15 CNTL. The first parameter gives the type of the matrix: 0 for unsymmetric matrix, 1 for symmetric positive matrix and 2 for general symmetric. The second parameter defined if the host processor work during the factorization and solves steps : PAR=1 host processor working and PAR=0 host processor not working.
The parameter INCTL and CNTL is the control parameter of MUMPS. The vectors ICNTL and CNTL in MUMPS becomes with index 1 like vector in fortran. A short description of all parameters of ICNTL and CNTL is given in ffmumps\_fileparam.txt $\codered$. For more details see the users' guide \cite{mumpsuserguide} $\codered$.

We describe now some elements of the main parameters of ICNTL for MUMPS.

* __Input matrix parameter__
	The input matrix is controlled by parameters ICNTL(5) and ICNTL(18). The matrix format (resp. matrix pattern and matrix entries) are controlled by INCTL(5) (resp. INCTL(18)). The different values of ICNTL(5) are 0 for assembled format and 1 for element format. In the current release of Freefem++, we consider that FE matrix or matrix is storage in assembled format. Therefore, INCTL(5) is treated as 0 value. The main option for ICNTL(18): INCLTL(18)=0 centrally on the host processor, ICNTL(18)=3 distributed the input matrix pattern and the entries (recommended option for distributed matrix by developer of MUMPS). For other values of ICNTL(18) see the user's guide of MUMPS. These values can be used also in FreeFem++.

	The default option implemented in FreeFem++ are ICNTL(5)=0 and ICNTL(18)=0.

* __Preprocessing parameter__
	The preprocessed matrix $A_{p}$ that will be effectively factored is defined by
	$$
	A_{p} = P \: D_r \: A \: Q_c \ D_c P^t
	$$
	where $P$ is the permutation matrix, $Q_c$ is the column permutation, $D_r$ and $D_c$ are diagonal matrix for respectively row and column scaling. The ordering strategy to obtain $P$ is controlled by parameter ICNTL(7). The permutation of zero free diagonal $Q_c$ is controlled by parameter ICNTL(6). The row and column scaling is controlled by parameter ICNTL(18). These option are connected and also strongly related with ICNTL(12) (see documentation of mumps for more details \cite{mumpsuserguide} $\codered$). The parameters `:::freefem permr`, `:::freefem scaler`, and `:::freefem scalec` in FreeFem++ allow to give permutation matrix($P$), row scaling ($D_r$) and column scaling ($D_c$) of the user respectively.

__Calling MUMPS in FreeFem++__

To call MUMPS in FreeFem++, we need to load the dynamic library MUMPS\_freefem.dylib (MacOSX), MUMPS\_freefem.so (Unix) or MUMPS\_freefem.dll (Windows) $\codered$.
This is done in typing `:::freefem load "MUMPS\_freefem"` in the .edp file. We give now the two methods to give the option of MUMPS solver in FreeFem++.

* __Solver parameters is defined in .edp file:__
	In this method, we need to give the parameters `:::freefem lparams` and `:::freefem dparams`. These parameters are defined for MUMPS by :

	`:::freefem lparams[0] = SYM`,<br>
	`:::freefem lparams[1] = PAR`,<br>
	$\forall i$ = 1,...,40, `:::freefem lparams[i+1] = ICNTL(i)`<br>
	$\forall i$ = 1,...,15, `:::freefem dparams[i-1] = CNTL(i)`

* __Reading solver parameters on a file:__

	The structure of data file for MUMPS in FreeFem++ is : first line parameter SYM and second line parameter PAR and in the following line the different value of vectors ICNTL and CNTL. An example of this parameter file is given in `:::freefem ffmumpsfileparam.txt`.

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

__Example__

A simple example of calling MUMPS in FreeFem++ with this two methods is given in the file testsolver_MUMPS.edp $\codered$ in the directory examples++-mpi.

### SuperLU distributed solver

The package SuperLU_DIST \cite{slu2,slu1} $\codered$ solves linear systems using LU factorization. It is a free scientific library under BSD license. The web site of this project is [http://crd.lbl.gov/~xiaoye/SuperLU](http://crd.lbl.gov/~xiaoye/SuperLU). This library provides functions to handle square or rectangular matrix in real and complex arithmetics. The method implemented in SuperLU_DIST is a supernodal method \cite{slu1} $\codered$. New release of this package includes a parallel symbolic factorization \cite{slu2} $\codered$. This scientific library is written in C and MPI for communications.

__Installation of SuperLU_DIST:__

To use SuperLU_DIST in FreeFem++, you have to install SuperLU_DIST package. We need MPI and ParMetis library to do this compilation. An installation procedure to obtain this package is given in the file README\_COMPILE in the directory src/solver/ $\codered$ of the FreeFem++ package.

__Creating Library of SuperLU_DIST interface for FreeFem++:__

The FreeFem++ interface to SuperLU_DIST for real (resp. complex) arithmetics is given in file real_SuperLU_DIST_FreeFem.cpp (resp. complex_SuperLU_DIST_FreeFem.cpp). These files are in the directory src/solver/. These interfaces are compatible with the release 3.2.1 of SuperLU_DIST. To use SuperLU_DIST in FreeFem++, we need libraries corresponding to these interfaces. A description to obtain these libraries is given in the file README_COMPILE in the directory src/solver $\codered$ of FreeFem++. We recall here the procedure. Go to the directory src/solver in FreeFem++ package. Edit the file makefile-sparsesolver.inc in your system : comment Section 1, comment line corresponding to libraries BLAS, Metis, ParMetis in Section 2 and comment in Section 3 the paragraph corresponding to SuperLU_DIST solver. And just type `:::bash make rsludist` (resp. `:::bash make csludist`) in the terminal to obtain the dynamic library of interface for real (resp. complex) arithmetics.

Now we give a short description of SuperLU_DIST parameters before describing the method to call SuperLU_DIST in FreeFem++.

__SuperLU_DIST parameters:__

We describe now some parameters of SuperLU_DIST. The SuperLU_DIST library use a 2D-logical process group. This process grid is specified by $nprow$ (process row) and $npcol$ (process column) such that $N_{p} = nprow \: npcol$ where $N_{p}$ is the number of all process allocated for SuperLU_DIST.

The input matrix parameters is controlled by "matrix= " in sparams for internal parameter or in the third line of parameters file. The different value are

* `:::freefem matrix = assembled` global matrix are available on all process
* `:::freefem matrix = distributedglobal` The global matrix is distributed among all the process
* `:::freefem matrix = distributed` The input matrix is distributed (not yet implemented)

The option arguments of SuperLU_DIST are described in the section Users-callable routine of \cite{sluuserguide} $\codered$. The parameter Fact and TRANS are specified in FreeFem++ interfaces to SuperLU_DIST during the different steps. For this reason, the value given by the user for this option is not considered.

The factorization LU is calculated in SuperLU_DIST on the matrix $A_p$.
$$
A_{p} = P_{c} \: P_r \: D_r \: A \: D_{c} \: P_{c}^{t}
$$
where $P_c$ and $P_r$ is the row and column permutation matrix respectively, $D_r$ and $D_c$ are diagonal matrix for respectively row and column scaling.
The option argument RowPerm (resp. ColPerm) control the row (resp. column) permutation matrix. $D_r$ and $D_c$ is controlled by the parameter DiagScale.
The parameter `:::freefem permr`, `:::freefem permc`, `:::freefem scaler`, and `:::freefem scalec` in FreeFem++ is provided to give row permutation, column permutation, row scaling and column scaling of the user respectively.
The other parameters for LU factorization are ParSymFact and ReplaceTinyPivot. The parallel symbolic factorization works only on a power of two processes and need the ParMetis ordering \cite{parmetis} $\codered$. The default option argument of SuperLU_DIST are given in the file ffsuperlu_dist_fileparam.txt.

__Calling SuperLU_DIST in FreeFem++__

To call SuperLU_DIST in FreeFem++, we need to load the library dynamic correspond to interface.
This done by the following line `:::freefem load "real_superlu _DIST_FreeFem"` (resp. `:::freefem load "complex_superlu_DIST_FreeFem"`) for real (resp. complex) arithmetics in the file .edp.

__Solver parameters is defined in .edp file:__

To call SuperLU_DIST with internal parameter, we used the parameters sparams. The value of parameters of SuperLU_DIST in sparams is defined by :

* nprow=1,
* npcol=1,
* matrix= distributedgloba,
* Fact= DOFACT,
* Equil=NO,
* ParSymbFact=NO,
* ColPerm= MMD\_AT\_PLUS\_A,
* RowPerm= LargeDiag,
* DiagPivotThresh=1.0,
* IterRefine=DOUBLE,
* Trans=NOTRANS,
* ReplaceTinyPivot=NO,
* SolveInitialized=NO,
* PrintStat=NO,
* DiagScale=NOEQUIL

This value correspond to the parameter in the file ffsuperlu_dist_fileparam.txt. $\codered$ If one parameter is not specified by the user, we take the default value of SuperLU_DIST.

__Reading solver parameters on a file:__
The structure of data file for SuperLU_DIST in FreeFem++ is given in the file ffsuperlu_dist_fileparam.txt (default value of the FreeFem++ interface).

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

__Example__

A simple example of calling SuperLU_DIST in FreeFem++ with this two methods is given in the file testsolver_superLU_DIST.edp in the directory examples++-mpi. $\codered$

### Pastix solver

Pastix (Parallel Sparse matrix package) is a free scientific library under CECILL-C license. This package solves sparse linear system with a direct and block ILU(k) iterative methods. This solver can be applied to a real or complex matrix with a symmetric pattern \cite{pastix} $\codered$.

__Installation of Pastix:__

To used Pastix in FreeFem++, you have to install pastix package in first. To compile this package, we need a fortran 90 compiler, scotch \cite{scotch} $\codered$ or Metis \cite{metis} $\codered$ ordering library and MPI. An installation procedure to obtain this package is given in the file .src/solver/ README_COMPILE in the section pastix of the FreeFem++ package.

__Creating Library of pastix interface for FreeFem++:__

The FreeFem++ interface to pastix is given in file real_pastix_FreeFem.cpp (resp. complex_pastix_FreeFem.cpp) for real (resp.complex) arithmetics. This interface is compatible with the release 2200 of pastix and is designed for a global matrix. We have also implemented interface for distributed matrices. To use pastix in FreeFem++, we need the library corresponding to this interface. A description to obtain this library is given in the file README_COMPILE in the directory src/solver of FreeFem++. We recall here the procedure. Go to the directory src/solver in FreeFem++ package. Edit the file makefile-sparsesolver.inc to yours system : comment Section 1, comment line corresponding to libraries BLAS, METIS and SCOTCH in Section 2 and comment in Section 3 the paragraph corresponding to pastix solver. And just type `:::bash make rpastix` (resp. `:::bash make cpastix`) in the terminal to obtain the dynamic library of interface for real (resp. complex) arithmetics.

Now we give a short description of pastix parameters before describing the method to call pastix in FreeFem++.

__Pastix parameters: __

The input `:::freefem matrix` parameter of FreeFem++ depend on pastix interface. `:::freefem matrix = assembled` for non distributed matrix. It is the same parameter for SuperLU_DIST. There are four parameters in Pastix : iparm, dparm, perm and invp. These parameters are respectively the integer parameters (vector of size 64), real parameters (vector of size 64), permutation matrix and inverse permutation matrix respectively. iparm and dparm vectors are described in \cite{pastixrefcard} $\codered$.
The parameters `:::freefem permr` and `:::freefem permc` in FreeFem++ are provided to give permutation matrix and inverse permutation matrix of the user respectively.

__Solver parameters defined in .edp file:__

To call Pastix in FreeFem++ in this case, we need to specify the parameters __lparams__ and __dparams__. These parameters are defined by :
$\codered$

$\forall i$ = 0,... ,63, `:::freefem lparams[i] = iparm[i]`.

$\forall i$ = 0,... ,63, `:::freefem dparams[i] = dparm[i]`.


__Reading solver parameters on a file:__

The structure of data file for pastix parameters in FreeFem++ is : first line structure parameters of the matrix and in the following line the value of vectors iparm and dparm in this order.

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

An example of this file parameter is given in ffpastix_iparm_dparm.txt with a description of these parameters. This file is obtained with the example file iparm.txt and dparm.txt including in the pastix package.

If no solver parameter is given, we use the default option of pastix solver.

__Example:__
A simple example of calling pastix in FreeFem++ with this two methods is given in the file testsolver_pastix.edp in the directory examples++-mpi. $\codered$

In Table 11.3 \ref{recap.direct.solveur} $\codered$, we recall the different matrix considering in the different direct solvers.

<table>
	<thead>
		<tr>
			<th colspan="7">Table 11.3 : Type of matrix used by the different direct sparse solver</th>
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

Concerning __iterative solvers__, we have chosen _pARMS_ \cite{spARMS} $\codered$, _HIPS_ \cite{HIPS} $\codered$ and _Hypre_ \cite{HYPRE} $\codered$.
Each software implements a different type of parallel preconditioner.
So, _pARMS_ implements algebraic domain decomposition preconditioner type such as additive Schwartz \cite{CAI89} $\codered$ and interface method \cite{LISAAD} $\codered$; while HIPS implement hierarchical incomplete factorization \cite{A:LaBRI::sisc06} $\codered$ and finally HYPRE implements multilevel preconditioner are AMG(Algebraic MultiGrid) \cite{A:AMG::sisc00} $\codered$ and parallel approximated inverse \cite{A:PASA::sisc00} $\codered$.

To use one of these programs in FreeFem++, you have to install it independently
of FreeFem++. It is also necessary to install the MPI communication library which is essential
for communication between the processors and, in some cases, software partitioning graphs like
METIS \cite{metis} $\codered$ or Scotch \cite{scotch} $\codered$.

All this preconditioners are used with Krylov subspace methods accelerators.
Krylov subspace methods are iterative methods which consist in finding a solution $x$ of linear system $Ax=b$ inside the affine space $x_0+K_m$ by imposing that $b-Ax \bot \mathcal{L}_m$, where $K_m$ is Krylov subspace of dimension $m$ defined by $K_m=\{r_0, Ar_0, A^2r_0,...,A^{m-1}r_0\}$ and $\mathcal{L}_m$ is another subspace of dimension $m$ which depends on type of Krylov subspace. For example in GMRES, $\mathcal{L}_m=AK_m$.

We realized an interface which is easy to use, so that the call of these different softwares in FreeFem++ is done in the same way. You just have to load the solver and then specify the parameters to apply to the specific solvers. In the rest of this chapter, when we talk about Krylov subspace methods we mean one among GMRES, CG and BICGSTAB.

### pARMS solver

_pARMS_ (parallel Algebraic Multilevel Solver) is a software developed by Youssef Saad and al at University of Minnesota \cite{spARMS} $\codered$.
This software is specialized in the resolution of large sparse non symmetric linear systems of equation. Solvers developed in pARMS is the Krylov subspace type.
It consists of variants of GMRES like FGMRES (Flexible GMRES) , DGMRES (Deflated GMRES) \cite{SAAD03} $\codered$ and BICGSTAB. pARMS also implements parallel preconditioner like RAS (Restricted Additive Schwarz)\cite{CAI89} $\codered$ and Schur Complement type preconditioner \cite{LISAAD} $\codered$.

All these parallel preconditioners are based on the principle of domain decomposition. Thus, the matrix $A$ is partitioned into sub matrices $A_i$($i=1,...,p$) where p represents the number of partitions one needs. The union of $A_i$ forms the original matrix. The solution of the overall system is obtained by solving the local systems on $A_i$ (see \cite{Smith96} $\codered$).
Therefore, a distinction is made between iterations on $A$ and the local iterations on $A_i$.
To solve the local problem on $A_i$ there are several preconditioners as __ilut__ (Incomplete LU with threshold), __iluk__ (Incomplete LU with level of fill in) and __ARMS__ (Algebraic Recursive Multilevel Solver). But to use pAMRS in FreeFem++ you have first to install pAMRS.

__Installation of pARMS__

To install pARMS, you must first download the pARMS package at \cite{spARMS} $\codered$. Once the download is complete, you must unpack package pARMS and follow the installation procedure described in file README to create the library __libparms.a__.

__Using pARMS as interface to FreeFem++__
Before calling pARMS solver inside FreeFem++, you must compile file parms\_FreeFem.cpp to create a dynamic library parms\_FreeFem.so. To do this, move to the directory src/solver of FreeFem++, edit the file makefile\-parms.inc to specify the following variables: $\codered$

```freefem
PARMS_DIR : // Directory of pARMS
PARMS_INCLUDE : // Directory for header of pARMS
METIS : // METIS directory
METIS_LIB : // METIS librairy
MPI : // MPI directory
MPI_INCLUDE : // MPI headers
FREEFEM : // FreeFem++ directory
FREEFEM_INCLUDE : // FreeFem++ header for sparse linear solver
LIBBLAS : // Blas library
```

After that, in the command line type `:::bash make parms` to create parms\_FreeFem.so.

As usual in FreeFem++, we will show by examples how to call pARMS in FreeFem++.
There are three ways of doing this:

__Example 1: Default parameters__
This example comes from user guide of FreeFem++ \cite{ufreefem} $\codered$ at page 12.

```freefem
load parms_freefem // Tell FreeFem that you will use pARMS
border C(t=0,2*pi){x=cos(t); y=sin(t);label=1;}
mesh Th = buildmesh (C(50));
fespace Vh(Th,P2);
Vh u,v;
func f= x*y;
problem Poisson(u,v,solver=sparsesolver) = // Bilinear part will use
int2d(Th)(dx(u)*dx(v) + dy(u)*dy(v)) // A sparse solver, in this case pARMS
- int2d(Th)( f*v) // Right hand side
+ on(1,u=0); // Dirichlet boundary condition

real cpu=clock();
Poisson; // SOLVE THE PDE
plot(u);
cout << " CPU time = " << clock()-cpu << endl;
```

In line 1 of example \ref{exm:first} $\codered$ we load in memory the pARMS dynamic library
with interface FreeFem++. After this, in line 7 we specify that the bilinear form will be solved by the
last sparse linear solver load in memory which, in this case, is pARMS.

The parameter used in pARMS in this case is the default one since the user does not have to provide any parameter.

Here are some default parameters:

* solver=FGMRES,
* Krylov dimension=30,
* Maximum of Krylov=1000,
* Tolerance for convergence=$1e-08$.(see book of Saad \cite{SAAD03} $\codered$ to
understand all this parameters.)
* preconditionner=Restricted Additif Schwarz \cite{CAI89} $\codered$,
* Inner Krylov dimension=5,
* Maximum of inner Krylov dimension=5,
* Inner preconditionner=ILUK.

To specify the parameters to apply to the solver, the user can either give an
integer vector for __integer parameters__ and real vectors for __real
parameters__ or provide a __file__ which contains those parameters.

__Example 2: User specifies parameters inside two vectors__

Lets us consider Navier Stokes example \ref{exm:segond} $\codered$ . In this example we
solve linear systems coming from discretization of Navier Stokes equation with pARMS. Parameters of solver is specified by user.

__Example 11.4 Stokes.edp__

```freefem
include "manual.edp"
include "includes.edp";
include "mesh_with_cylinder.edp";
include "bc_poiseuille_in_square.edp";
include "fe_functions.edp";
load parms_FreeFem
int[int] iparm(16); real[int] dparm(6);
int ,ii;
for(ii=0;ii<16;ii++){iparm[ii]=-1;} for(ii=0;ii<6;ii++) dparm[ii]=-1.0;
fespace Vh(Th,[P2,P2,P1]);
iparm[0]=0;
varf Stokes ([u,v,p],[ush,vsh,psh],\textbf{solver=sparsesolver}) =
 int2d(Th)( nu*( dx(u)*dx(ush) + dy(u)*dy(ush) + dx(v)*dx(vsh) + dy(v)*dy(vsh) )
 - p*psh*(1.e-6) 			 // p epsilon
 - p*(dx(ush)+dy(vsh)) //+ dx(p)*ush + dy(p)*vsh
 - (dx(u)+dy(v))*psh 	 // psh div(u)
 )
 + on(cylinder,infwall,supwall,u=0.,v=0.)+on(inlet,u=uc,v=0); // Bdy conditions
matrix AA=Stokes(VVh,VVh);
set(AA,solver=sparsesolver,lparams=iparm,dparams=dparm); //Set pARMS as linear solver
real[int] bb= Stokes(0,VVh); real[int] sol(AA.n);
sol= AA^-1 * bb;
```

We need two vectors to specify the parameters of the linear solver. In line 1 of example \ref{exm:segond} $\codered$ we have declared these vectors(`:::freefem int[int] iparm(16); real[int] dparm(6);`). In line 3 we have initialized these vectors by negative values. We do this because all parameters values in pARMS are positive and if you do not change the negative values of one entry of this vector, the default value will be set. In tables (table 11.4 \ref{lpparm} $\codered$ and 11.5 \ref{pardoubleparm} $\codered$), we have the meaning of differents entries of these vectors.

<table>
	<thead>
		<tr>
			<th colspan="2">Table 11.4 : Meaning of lparams corresponding variables for example
			\ref{exm:segond} $\codered$</th>
		</tr>
	</thead>
	<tbody>
		<tr>
			<td style="font-weight: bold">Entries of iparm</td>
			<td style="font-weight: bold">Significations of each entries</td>
		</tr>
	    <tr>
	        <td>iparm[0]</td>
	        <td>Krylov subspace methods.<br>Differents valuesfor this parameters are specify on table \ref{kryparms} $\codered$</td>
	    </tr>
	    <tr>
	        <td>iparm[1]</td>
	        <td>Preconditionner.<br>Differentspreconditionners for this parameters are specify on table \ref{precond} $\codered$</td>
	    </tr>
	    <tr>
	        <td>iparm[2]</td>
	        <td>Krylov subspace dimension in outer iteration: default value 30</td>
	    </tr>
	    <tr>
	        <td>iparm[3]</td>
	        <td>Maximum of iterations in outer iteration: default value 1000</td>
	    </tr>
	    <tr>
	        <td>iparm[4]</td>
	        <td>Number of level in arms when used.</td>
	    </tr>
	    <tr>
	        <td>iparm[5]</td>
	        <td>Krylov subspace dimension in inner iteration: default value 3</td>
	    </tr>
	    <tr>
	        <td>iparm[6]</td>
	        <td>Maximum of iterations in inner iteration: default value 3</td>
	    </tr>
	    <tr>
	        <td>iparm[7]</td>
	        <td>Symmetric(=1 for symmetric) or unsymmetric matrix:<br>default value 0(unsymmetric matrix)</td>
	    </tr>
	    <tr>
	        <td>iparm[8]</td>
	        <td>Overlap size between different subdomain: default value 0(no overlap)</td>
	    </tr>
	    <tr>
	        <td>iparm[9]</td>
	        <td>Scale the input matrix or not: Default value 1 (Matrix should bescale)</td>
	    </tr>
	    <tr>
	        <td>iparm[10]</td>
	        <td>Block size in arms when used: default value 20</td>
	    </tr>
	    <tr>
	        <td>iparm[11]</td>
	        <td>lfil0 (ilut, iluk, and arms) : default value 20</td>
	    </tr>
	    <tr>
	        <td>iparm[12]</td>
	        <td>lfil for Schur complement const : default value 20</td>
	    </tr>
	    <tr>
	        <td>iparm[13]</td>
	        <td>lfil for Schur complement const : default value 20</td>
	    </tr>
	    <tr>
	        <td>iparm[14]</td>
	        <td>Multicoloring or not in ILU when used : default value 1</td>
	    </tr>
	    <tr>
	        <td>iparm[15]</td>
	        <td>Inner iteration : default value 0</td>
	    </tr>
	    <tr>
	        <td>iparm[16]</td>
	        <td>Print message when solving:default 0 (no messageprint).<br>0: no message is print,<br>1: Convergence informations like number of iteration and residual ,<br>2: Timing for a different step like preconditioner<br>3 : Print all informations.</td>
	    </tr>
	</tbody>
</table>

<table>
	<thead>
		<tr>
			<th colspan="2">Table 11.5 : Significations of dparams corresponding variables for example \ref{exm:segond} $\codered$</th>
		</tr>
	</thead>
	<tbody>
	    <tr>
	        <td style="font-weight: bold">Entries of dparm</td>
	        <td style="font-weight: bold">Significations of each entries</td>
	    </tr>
	    <tr>
	        <td>dparm[0]</td>
	        <td>precision for outer iteration : default value 1e-08</td>
	    </tr>
	    <tr>
	        <td>dparm[1]</td>
	        <td>precision for inner iteration: default value 1e-2</td>
	    </tr>
	    <tr>
	        <td>dparm[2]</td>
	        <td>tolerance used for diagonal domain: : default value 0.1</td>
	    </tr>
	    <tr>
	        <td>dparm[3]</td>
	        <td>drop tolerance droptol0 (ilut, iluk, and arms) : default value 1e-2</td>
	    </tr>
	    <tr>
	        <td>dparm[4]</td>
	        <td>droptol for Schur complement const: default value 1e-2</td>
	    </tr>
	    <tr>
	        <td>dparm[5]</td>
	        <td>droptol for Schur complement const: default value 1e-2</td>
	    </tr>
	</tbody>
</table>

<table>
	<thead>
		<tr>
			<th colspan="2">Table 11.6 : Krylov Solvers in pARMS</th>
		</tr>
	</thead>
	<tbody>
	    <tr>
	        <td style="font-weight: bold">Values of iparm[0]</td>
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
			<th colspan="2">Table 11.7 : Preconditionners in pARMS</th>
		</tr>
	</thead>
	<tbody>
	    <tr>
	        <td style="font-weight: bold">Values of iparm[1]</td>
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

We run example 11.4 \ref{exm:segond} $\codered$ on cluster paradent of Grid5000 and report results in table 11.8 \ref{parmResult} $\codered$.

<table>
	<thead>
		<tr>
			<th colspan="5">Table 11.8 : Convergence and time for solving linear system from example
			\ref{exm:segond} 11.4 $\codered$</th>
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
	        <td align="center">nit</td>
	        <td align="center">time</td>
	        <td align="center">nit</td>
	        <td align="center">time</td>
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
			<th colspan="2">Table 11.9 : Legend of table 11.8</th>
		</tr>
	</thead>
	<tbody>
    <tr>
        <td>n</td>
        <td>matrix size</td>
    </tr>
    <tr>
        <td>nnz</td>
        <td>number of non null entries inside matrix</td>
    </tr>
    <tr>
        <td>nit</td>
        <td>number of iteration for convergence</td>
    </tr>
    <tr>
        <td>time</td>
        <td>Time for convergence</td>
    </tr>
    <tr>
        <td>Te</td>
        <td>Time for constructing finite element matrix</td>
    </tr>
    <tr>
        <td>np</td>
        <td>number of processor</td>
    </tr>
</table>

In this example, we fix the matrix size (in term of finite element, we fix the mesh) and increase the number of processors used to solve the linear system. We saw that, when the number of processors increases, the time for solving the linear equation decreases, even if the number of iteration increases. This proves that, using pARMS as solver of linear systems coming from discretization of partial differential equation in FreeFem++ can decrease drastically the total time of simulation.
