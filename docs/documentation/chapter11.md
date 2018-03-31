# Parallel sparse solvers

Parallel sparse solvers use several processors to solve linear systems of equation. Like sequential,
parallel linear solvers can be direct or iterative. In FreeFem++ both are available.

## Using parallel sparse solvers in FreeFem++

We recall that the `:::freefem solver` parameters are defined in the following commands: `:::freefem solve`, `:::freefem problem`, `:::freefem set` (setting parameter of a matrix) and in the construction of the matrix corresponding to a bilinear form. In these commands, the parameter `:::freefem solver` must be set to `:::freefem sparsesolver` for parallel sparse solver. We have added specify parameters to these command lines for parallel sparse solvers. These are

* `:::freefem lparams:` vector of integer parameters (l is for the c++ type long)
* `:::freefem dparams:` vector of real parameters
* `:::freefem sparams:` string parameters
* `:::freefem datafilename:` name of the file which contains solver parameters

The following four parameters are only for direct solvers and are vectors. These parameters allow the user to preprocess the matrix (see the section on sparse direct solver above for more information).

* `:::freefem permr:` row permutation (integer vector)
* `:::freefem permc:` column permutation or inverse row permutation (integer vector)
* `:::freefem scaler:` row scaling (real vector)
* `:::freefem scalec:` column scaling (real vector)


There are two possibilities to control solver parameters. The first method defines parameters with `:::freefem lparams`, `:::freefem dparams` and `:::freefem sparams` in .edp file. The second one reads the solver parameters from a data file. The name of this file is specified by `:::freefem datafilename`. If `:::freefem lparams`, `:::freefem dparams`, `:::freefem sparams` or `:::freefem datafilename` is not provided by the user, the solver's default value is used.

To use parallel solver in FreeFem++, we need to load the dynamic library corresponding to this solver. For example to use MUMPS solver as parallel solver in FreeFem, write in the .edp file __load "MUMPS\_FreeFem"__.

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
A description to obtain this library is given in the file README\_COMPILE in the directory src/solver $\codered$ of FreeFem++. We recall here the procedure. Go to the directory src/solver in FreeFem++ package. Edit the file makefile-sparsesolver.inc to yours system: comment Section 1, comment line corresponding to libraries BLAS, BLACS, ScaLAPACK, Metis, scotch in Section 2 and comment in Section 3 $\codered$ the paragraph corresponding to MUMPS solver. And then type __make mumps__ in a terminal window.

Now we give a short description of MUMPS parameters before describing the method to call MUMPS in FreeFem++.

__MUMPS parameters:__

There are four input parameters in MUMPS (see \cite{mumpsuserguide} $\codered$). Two integers SYM and PAR, a vector of integer of size 40 INCTL and a vector of real of size 15 CNTL. The first parameter gives the type of the matrix: 0 for unsymmetric matrix, 1 for symmetric positive matrix and 2 for general symmetric. The second parameter defined if the host processor work during the factorization and solves steps : PAR=1 host processor working and PAR=0 host processor not working.
The parameter INCTL and CNTL is the control parameter of MUMPS. The vectors ICNTL and CNTL in MUMPS becomes with index 1 like vector in fortran. A short description of all parameters of ICNTL and CNTL is given in ffmumps\_fileparam.txt $\codered$. For more details see the users' guide \cite{mumpsuserguide} $\codered$.

We describe now some elements of the main parameters of ICNTL for MUMPS.

* __Input matrix parameter__
	The input matrix is controlled by parameters ICNTL(5) and ICNTL(18). The matrix format (resp. matrix pattern and matrix entries) are controlled by INCTL(5) (resp. INCTL(18)). The different values of ICNTL(5) are 0 for assembled format and 1 for element format. In the current release of Freefem$++$, we consider that FE matrix or matrix is storage in assembled format. Therefore, INCTL(5) is treated as 0 value. The main option for ICNTL(18): INCLTL(18)=0 centrally on the host processor, ICNTL(18)=3 distributed the input matrix pattern and the entries (recommended option for distributed matrix by developer of MUMPS). For other values of ICNTL(18) see the user's guide of MUMPS. These values can be used also in FreeFem++.

	The default option implemented in FreeFem++ are ICNTL(5)=0 and ICNTL(18)=0.

* __Preprocessing parameter__
	The preprocessed matrix $A_{p}$ that will be effectively factored is defined by
	$$
	A_{p} = P \: D_r \: A \: Q_c \ D_c P^t
	$$
	where $P$ is the permutation matrix, $Q_c$ is the column permutation, $D_r$ and $D_c$ are diagonal matrix for respectively row and column scaling. The ordering strategy to obtain $P$ is controlled by parameter ICNTL(7). The permutation of zero free diagonal $Q_c$ is controlled by parameter ICNTL(6). The row and column scaling is controlled by parameter ICNTL(18). These option are connected and also strongly related with ICNTL(12) (see documentation of mumps for more details \cite{mumpsuserguide} $\codered$). The parameters permr, scaler, and scalec in FreeFem++ allow to give permutation matrix($P$), row scaling ($D_r$) and column scaling ($D_c$) of the user respectively.

__Calling MUMPS in FreeFem++__

To call MUMPS in FreeFem++, we need to load the dynamic library MUMPS\_freefem.dylib (MacOSX), MUMPS\_freefem.so (Unix) or MUMPS\_freefem.dll (Windows) $\codered$.
This is done in typing load "MUMPS\_freefem" in the .edp file. We give now the two methods to give the option of MUMPS solver in FreeFem++.

* __Solver parameters is defined in .edp file:__
	In this method, we need to give the parameters `:::freefem lparams` and `:::freefem dparams`. These parameters are defined for MUMPS by :
	$\codered$
	\begin{center}
	\begin{tabular}{ll}
	&lparams[0] = SYM, \\
	&lparams[1] = PAR, \\
	$\forall i=1,\ldots,40, \quad$ & lparams[$i$+1] = ICNTL($i$).\\
	\\
	$\forall i=1,\ldots,15, \quad$ & dparams[$i-1$] = CNTL($i$).
	\end{tabular}
	\end{center}

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

The package SuperLU_DIST \cite{slu2,slu1} $\codered$ solves linear systems using LU factorization. It is a free scientific library under BSD license. The web site of this project is http://crd.lbl.gov/~xiaoye/SuperLU. This library provides functions to handle square or rectangular matrix in real and complex arithmetics. The method implemented in SuperLU_DIST is a supernodal method \cite{slu1} $\codered$. New release of this package includes a parallel symbolic factorization \cite{slu2} $\codered$. This scientific library is written in C and MPI for communications.

__Installation of SuperLU_DIST:__

To use SuperLU_DIST in FreeFem++, you have to install SuperLU_DIST package. We need MPI and ParMetis library to do this compilation. An installation procedure to obtain this package is given in the file README\_COMPILE in the directory src/solver/ of the FreeFem++ package.

__Creating Library of SuperLU_DIST interface for FreeFem++:__

The FreeFem++ interface to SuperLU_DIST for real (resp. complex) arithmetics is given in file real_SuperLU_DIST_FreeFem.cpp (resp. complex_SuperLU_DIST_FreeFem.cpp). These files are in the directory src/solver/. These interfaces are compatible with the release 3.2.1 of SuperLU_DIST. To use SuperLU_DIST in FreeFem++, we need libraries corresponding to these interfaces. A description to obtain these libraries is given in the file README_COMPILE in the directory src/solver of FreeFem++. We recall here the procedure. Go to the directory src/solver in FreeFem++ package. Edit the file makefile-sparsesolver.inc in your system : comment Section 1, comment line corresponding to libraries BLAS, Metis, ParMetis in Section 2 and comment in Section 3 the paragraph corresponding to SuperLU_DIST solver. And just type __make rsludist__ (resp. __make csludist__) in the terminal to obtain the dynamic library of interface for real (resp. complex) arithmetics.

Now we give a short description of SuperLU_DIST parameters before describing the method to call SuperLU_DIST in FreeFem++.

__SuperLU_DIST parameters:__

We describe now some parameters of SuperLU_DIST. The SuperLU_DIST library use a 2D-logical process group. This process grid is specifies by $nprow$ (process row) and $npcol$ (process column) such that $N_{p} = nprow \: npcol$ where $N_{p}$ is the number of all process allocated for SuperLU_DIST.

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
The parameter permr, permc, scaler, and scalec in FreeFem++ is provided to give row permutation, column permutation, row scaling and column scaling of the user respectively.
The other parameters for LU factorization are ParSymFact and ReplaceTinyPivot. The parallel symbolic factorization works only on a power of two processes and need the ParMetis ordering \cite{parmetis} $\codered$. The default option argument of SuperLU_DIST are given in the file ffsuperlu_dist_fileparam.txt.

__Calling SuperLU_DIST in FreeFem++__

To call SuperLU_DIST in FreeFem++, we need to load the library dynamic correspond to interface.
This done by the following line __load "real_superlu _DIST_FreeFem"__ (resp. __load "complex_superlu_DIST_FreeFem"__) for real (resp. complex) arithmetics in the file .edp.

__Solver parameters is defined in .edp file:__

To call SuperLU_DIST with internal parameter, we used the parameters sparams. The value of parameters of SuperLU_DIST in sparams is defined by
$\codered$
\begin{tabular}{ll}
sparams&="nprow=1, npcol=1, matrix= distributedgloba, Fact= DOFACT, Equil=NO, \\
	 & ParSymbFact=NO, ColPerm= MMD\_AT\_PLUS\_A, RowPerm= LargeDiag, \\
	 & DiagPivotThresh=1.0, IterRefine=DOUBLE, Trans=NOTRANS, \\
	 & ReplaceTinyPivot=NO, SolveInitialized=NO, PrintStat=NO, DiagScale=NOEQUIL "
\end{tabular}
This value correspond to the parameter in the file ffsuperlu_dist_fileparam.txt. If one parameter is not specify by the user, we take the default value of SuperLU_DIST.

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

A simple example of calling SuperLU_DIST in FreeFem++ with this two methods is given in the file testsolver_superLU_DIST.edp in the directory examples++-mpi.

### Pastix solver

Pastix (Parallel Sparse matrix package) is a free scientific library under CECILL-C license. This package solves sparse linear system with a direct and block ILU(k) iterative methods. This solver can be applied to a real or complex matrix with a symmetric pattern \cite{pastix} $\codered$.

__Installation of Pastix:__

To used Pastix in FreeFem++, you have to install pastix package in first. To compile this package, we need a fortran 90 compiler, scotch \cite{scotch} $\codered$ or Metis \cite{metis} $\codered$ ordering library and MPI. An installation procedure to obtain this package is given in the file .src/solver/ README_COMPILE in the section pastix of the FreeFem++ package.

__Creating Library of pastix interface for FreeFem++:__

The FreeFem++ interface to pastix is given in file real_pastix_FreeFem.cpp (resp. complex_pastix_FreeFem.cpp) for real (resp.complex) arithmetics. This interface is compatible with the release 2200 of pastix and is designed for a global matrix. We have also implemented interface for distributed matrices. To use pastix in FreeFem++, we need the library corresponding to this interface. A description to obtain this library is given in the file README_COMPILE in the directory src/solver of FreeFem++. We recall here the procedure. Go to the directory src/solver in FreeFem++ package. Edit the file makefile-sparsesolver.inc to yours system : comment Section 1, comment line corresponding to libraries BLAS, METIS and SCOTCH in Section 2 and comment in Section 3 the paragraph corresponding to pastix solver. And just type __make rpastix__ (resp. __make cpastix__) in the terminal to obtain the dynamic library of interface for real (resp. complex) arithmetics.

Now we give a short description of pastix parameters before describing the method to call pastix in FreeFem++.

__Pastix parameters: __

The input `:::freefem matrix` parameter of FreeFem++ depend on pastix interface. `:::freefem matrix = assembled` for non distributed matrix. It is the same parameter for SuperLU_DIST. There are four parameters in Pastix : iparm, dparm, perm and invp. These parameters are respectively the integer parameters (vector of size 64), real parameters (vector of size 64), permutation matrix and inverse permutation matrix respectively. iparm and dparm vectors are described in \cite{pastixrefcard} $\codered$.
The parameters permr and permc in FreeFem++ are provided to give permutation matrix and inverse permutation matrix of the user respectively.

__Solver parameters defined in .edp file:__

To call Pastix in FreeFem++ in this case, we need to specify the parameters __lparams__ and __dparams__. These parameters are defined by :
$\codered$
\begin{center}
\begin{tabular}{ll}
$\forall i=0,\ldots,63, \quad$ & lparams[$i$] = iparm[$i$].\\
\\
$\forall i=0,\ldots,63, \quad$ & dparams[$i$] = dparm[$i$].
\end{tabular}
\end{center}


__Reading solver parameters on a file:__

The structure of data file for pastix parameters in FreeFem++ is : first line structure parameters of the matrix and in the following line the value of vectors iparm and dparm in this order.

```freefem
assembled /* matrix input :: assembled, distributed global and distributed */
iparm[0]
iparm[1]
\ldots
\ldots
iparm[63]
dparm[0]
dparm[1]
\ldots
\ldots
dparm[63]
```

An example of this file parameter is given in ffpastix_iparm_dparm.txt with a description of these parameters. This file is obtained with the example file iparm.txt and dparm.txt including in the pastix package.

If no solver parameter is given, we use the default option of pastix solver.

__Example:__
A simple example of calling pastix in FreeFem++ with this two methods is given in the file testsolver_pastix.edp in the directory examples++-mpi.

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
It consists of variants of GMRES like FGMRES(Flexible GMRES) , DGMRES(Deflated GMRES) \cite{SAAD03} $\codered$ and BICGSTAB. pARMS also implements parallel preconditioner like RAS (Restricted Additive Schwarz)\cite{CAI89} $\codered$ and Schur Complement type preconditioner \cite{LISAAD} $\codered$.

All these parallel preconditioners are based on the principle of domain decomposition. Thus, the matrix $A$ is partitioned into sub matrices $A_i$($i=1,...,p$) where p represents the number of partitions one needs. The union of $A_i$ forms the original matrix. The solution of the overall system is obtained by solving the local systems on $A_i$ (see \cite{Smith96} $\codered$).
Therefore, a distinction is made between iterations on $A$ and the local iterations on $A_i$.
To solve the local problem on $A_i$ there are several preconditioners as __ilut__ (Incomplete LU with threshold), __iluk__ (Incomplete LU with level of fill in) and __ARMS__ (Algebraic Recursive Multilevel Solver). But to use pAMRS in FreeFem++ you have first to install pAMRS.

__Installation of pARMS__

To install pARMS, you must first download the pARMS package at \cite{spARMS} $\codered$. Once the download is complete, you must unpack package pARMS and follow the installation procedure described in file README to create the library __libparms.a__.

__Using pARMS as interface to FreeFem++__
Before calling pARMS solver inside FreeFem++, you must compile file $parms\_FreeFem.cpp$ to create a dynamic library $parms\_FreeFem.so$. To do this, move to the directory $src/solver$ of FreeFem++, edit the file $makefile\-parms.inc$ to specify the following variables:

\begin{tabular}{ll}
 %\begin{tabular*}
$PARMS\_DIR$ : & Directory of pARMS \\
\textbf{$PARMS\_INCLUDE$} : & Directory for header of pARMS \\
\textbf{$METIS$} : & METIS directory \\
\textbf{$METIS\_LIB$} : & METIS librairy \\
\textbf{$MPI$} : & MPI directory \\
\textbf{$MPI\_INCLUDE$} : & MPI headers \\
\textbf{$FREEFEM$} : & FreeFem++ directory \\
\textbf{$FREEFEM\_INCLUDE$} : & FreeFem++ header for sparse linear solver\\
\textbf{$LIBBLAS$} : & Blas library\\
\end{tabular}

After that, in the command line type `:::bash make parms` to create $parms\_FreeFem.so$.

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
problem Poisson(u,v,solver=sparsesolver) = // bilinear part will use
int2d(Th)(dx(u)*dx(v) + dy(u)*dy(v)) // a sparse solver, in this case pARMS
- int2d(Th)( f*v) // right hand side
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
			<th colspan="2">Table 11.4 : Meaning of __lparams__ corresponding variables for example
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
	        <td>Print message when solving:default 0(no messageprint).<br>0: no message is print,<br>1: Convergence informations like number of iteration and residual ,<br>2: Timing for a different step like preconditioner<br>3 : Print all informations.</td>
	    </tr>
	</tbody>
</table>

<table>
	<thead>
		<tr>
			<th colspan="2">Table 11.5 : Significations of __dparams__ corresponding variables for example \ref{exm:segond} $\codered$</th>
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

### Interfacing with HIPS

__HIPS__ (_Hierarchical Iterative Parallel Solver_) is a scientific library that provides an efficient parallel iterative solver for very large sparse linear systems. HIPS is available as free software under the CeCILL-C licence. The interface that we realized is compatible with release __1.2 beta.rc4__ $\codered$ of HIPS.

HIPS implements two solver classes which are the iteratives class (GMRES, PCG) and the Direct class. Concerning preconditionners, HIPS implements a type of multilevel ILU. For further informations on those preconditionners see \cite{A:LaBRI::sisc06,
A:LaBRI::HRR07} $\codered$.

__Installation of HIPS__

To install HIPS, first download the HIPS package at \cite{HIPS} $\codered$, unpack it and go to the HIPS source directory. The installation of HIPS is machine dependence. For example, to install HIPS on a linux cluster copy the file __$Makefile\_Inc\_Examples/makefile.inc.gnu$__ $\codered$ on the root directory of HIPS with the name __makefile.inc__. After this, edit __makefile.inc__ to set values of different variables and type `:::bash make all`.

__Using HIPS as the interface to FreeFem++__

Before calling the HIPS solver inside FreeFem++, you must compile file $hips\_FreeFem.cpp$ to create dynamic library $hips\_FreeFem.so$. To do this, move to the directory $src/solver$ of FreeFem++ and edit the file
$makefile.inc$ to specify the following variables: $\codered$

\begin{tabular}{ll}
 %\begin{tabular*}
\textbf{$HIPS\_DIR$} : & Directory of HIPS \\
\textbf{$HIPS\_INCLUDE$}: & -I\$($HIPS\_DIR$)/SRC/INCLUDE : Directory for HIPS
headers\\
\textbf{$LIB\_DIR$} : & -L\$($HIPS\_DIR$)/LIB : Librairies directory \\
\textbf{$LIBHIPSSEQUENTIAL$} : & \$($HIPS\_DIR$)/LIB/libhipssequential.a: HIPS
utilities library\\
\textbf{$LIBHIPS$} : & \$($HIPS\_DIR$)/LIB/libhips.a: HIPS library\\
\textbf{$FREEFEM$} : & FreeFem++ directory \\
\textbf{$FREEFEM\_INCLUDE$} : & FreeFem headers for sparse linear solver\\
\textbf{$METIS$} : & METIS directory \\
\textbf{$METIS\_LIB$} : & METIS library \\
\textbf{$MPI$} : & MPI directory \\
\textbf{$MPI\_INCLUDE$} : & MPI headers \\
\end{tabular}

After specifies all the variables, in the command line in the directory
$src/solver$ type \textbf{make hips} to create $hips\_FreeFem.so$. Like with pARMS, the calling of HIPS in FreeFem++ can be done in three different manners. We will present only one example where the user specifies the parameters through keywords `:::freefem lparams` and `:::freefem dparams`.

__Laplacian 3D solve with HIPS__

Let us consider the 3D Laplacian example inside FreeFem++ package where after discretization we want to solve the linear equation with Hips. Example \ref{hips:laplacian} $\codered$ is Laplacian3D using Hips as linear solver. We first load Hips solver at line 2. From line 4 to 15 we specify the parameters for the Hips solver and in line 46 of example \ref{hips:laplacian} $\codered$ we set these parameters in the linear solver.

In Table \ref{hipslabel} $\codered$ results of running example \ref{hips:laplacian} $\codered$ on Cluster Paradent of Grid5000 are reported. We can see in this running example the efficiency of parallelism.

__Example Laplacian3D.edp__

```freefem
load "msh3"
load "hips_FreeFem" //load library
int nn=10,iii;
int[int] iparm(14);
real[int] dparm(6);
for(iii=0;iii<14;iii++)iparm[iii]=-1;
for(iii=0;iii<6;iii++) dparm[iii]=-1;
iparm[0]=0; //use iterative solver
iparm[1]=1; //PCG as Krylov method
iparm[4]=0; //Matrix are symmetric
iparm[5]=1; //Pattern are also symmetric
iparm[9]=1; //Scale matrix
dparm[0]=1e-13;//Tolerance to convergence
dparm[1]=5e-4; //Threshold in ILUT
dparm[2]=5e-4; //Threshold for Schur preconditionner
mesh Th2=square(nn,nn);
 fespace Vh2(Th2,P2);
 Vh2 ux,uz,p2;
 int[int] rup=[0,2], rdown=[0,1], rmid=[1,1,2,1,3,1,4,1];
real zmin=0,zmax=1;
mesh3 Th=buildlayers(Th2,nn,
 zbound=[zmin,zmax],
 reffacemid=rmid,
 reffaceup = rup,
 reffacelow = rdown);
savemesh(Th,"copie.mesh");
mesh3 Th3("copie.mesh");
fespace Vh(Th,P2);
func ue = 2*x*x + 3*y*y + 4*z*z + 5*x*y+6*x*z+1;
func uex= 4*x+ 5*y+6*z;
func uey= 6*y + 5*x;
func uez= 8*z +6*x;
func f= -18. ;
Vh uhe = ue; //
cout << " uhe min: " << uhe[].min << " max:" << uhe[].max << endl;
Vh u,v;
macro Grad3(u) [dx(u),dy(u),dz(u)] // EOM
varf va(u,v)= int3d(Th)(Grad3(v)' *Grad3(u)) //') for emacs
 + int2d(Th,2)(u*v)
 - int3d(Th)(f*v)
 - int2d(Th,2) ( ue*v + (uex*N.x +uey*N.y +uez*N.z)*v )
 + on(1,u=ue);
real cpu=clock();
matrix Aa;
Aa=va(Vh,Vh);
varf l(unused,v)=int3d(Th)(f*v);
Vh F; F[]=va(0,Vh);
if(mpirank==0){
 cout << "Taille " << Aa.n << endl;
 cout << "Non zeros " << Aa.nbcoef << endl;
}
if(mpirank==0)
cout << "CPU TIME FOR FORMING MATRIX = " << clock()-cpu << endl;
set(Aa,solver=sparsesolver,dparams=dparm, lparams=iparm); //Set hips as linear solver
u[]=Aa^-1*F[];
```

\begin{table}
\begin{center}
\begin{tabular}{|c|c|c|}
\textbf{$n=4 \times 10^6$} & \textbf{$nnz=118 \times 10^6$} & \textbf{Te=221.34}
\\
np & nit & time \\
8 & 190 & 120.34 \\
16 & 189 & 61.08 \\
32 & 186 & 31.70 \\
64 & 183 & 23.44 \\
\end{tabular}
\end{center}
\caption{Iterations and Timing of solving linear system from example
\ref{hips:laplacian} $\codered$}

\end{table}
Legend of table \ref{hipslabel} $\codered$ are give in table \ref{legtableparm} $\codered$.

<table>
	<thead>
		<tr>
			<th colspan="2">Table 11.9 : Legend of table 11.8</th>
		</tr>
	</thead>
	<tbody>
	    <tr>
	        <td>$n=4 \times 10^6$</td>
	        <td>$nnz=118 \times 10^6$</td>
	        <td>Te=221.34</td>
	    </tr>
	    <tr>
	        <td class="align-center">np</td>
	        <td class="align-center">nit</td>
	        <td class="align-center">time</td>
	    </tr>
	    <tr>
	        <td class="align-center">8</td>
	        <td class="align-center">190</td>
	        <td class="align-center">120.34</td>
	    </tr>
	    <tr>
	        <td class="align-center">16</td>
	        <td class="align-center">189</td>
	        <td class="align-center">61.08</td>
	    </tr>
	    <tr>
	        <td class="align-center">32</td>
	        <td class="align-center">186</td>
	        <td class="align-center">31.70</td>
	    </tr>
	    <tr>
	        <td class="align-center">64</td>
	        <td class="align-center">183</td>
	        <td class="align-center">23.44</td>
	    </tr>
	</tbody>
</table>

\begin{table}[hbtp]
\begin{center}
\begin{tabular}{|l|l|}
\textbf{Entries of iparm } & \textbf{Significations of each entries} \\
\multirow{2}{*}{iparm[0] } & Strategy use for solving \\ &
 ( Iterative=0 or Hybrid=1 or Direct=2 ). Defaults values are : Iterative \\

\multirow{2}{*}{iparm[1] } & Krylov methods. \\ & If iparm[0]=0, give type of
Krylov methods: 0 for GMRES, 1 for PCG \\
iparm[2] & Maximum of iterations in outer iteration: default value 1000 \\

iparm[3] & Krylov subspace dimension in outer iteration: default value 40 \\


\multirow{2}{*}{iparm[4]} & Symmetric(=0 for symmetric) and 1 for unsymmetric
matrix: \\ & default value 1(unsymmetric matrix) \\
iparm[5] & Pattern of matrix are symmetric or not: default value 0 \\
iparm[6] & Partition type of input matrix: dafault value 0 \\
\multirow{2}{*}{iparm[7]} & Number of level that use the HIPS locally consistent
fill-in:\\ & Default value 2 \\
\multirow{2}{*}{iparm[8]} & Numbering in indices array will start at 0 or 1:\\
& Default value 0 \\
iparm[9] & Scale matrix. Default value 1 \\
\multirow{2}{*}{iparm[10]} & Reordering use inside subdomains for reducing
fill-in:\\ & Only use for iterative. Default value 1 \\
\multirow{2}{*}{iparm[11]} & Number of unknowns per node in the matrix non-zero
pattern graph: \\ & Default value 1 \\
\multirow{2}{*}{iparm[12]} & This value is used to set the number of time the
\\ & normalization is applied to the matrix: Default 2. \\
iparm[13] & Level of informations printed during solving: Default 5. \\
iparm[14] & HIPS\_DOMSIZE Subdomain size \\
\end{tabular}
\end{center}
\caption{Significations of \textbf{lparams} corresponding to HIPS interface }

\end{table}



\begin{table}[hbtp]
\begin{center}
\begin{tabular}{|l|l|}
%\begin{table}[hb]
dparm[0] & $HIPS\_PREC$: Relative residual norm: Default=1e-9 \\
\multirow{2}{*}\textbf{ dparm[1] }& $HIPS\_DROPTOL0$: Numerical threshold in
ILUT for interior domain \\
& (important : set 0.0 in HYBRID: Default=0.005)\\
\multirow{2}{*}{dparm[2]} & $HIPS\_DROPTOL1$ : Numerical threshold in ILUT
for\\
& Schur preconditioner: Default=0.005\\
\multirow{2}{*} \textbf{dparm[3] } & $HIPS\_DROPTOLE$ : Numerical threshold for
coupling between the \\
& interior level and Schur: Default 0.005\\

\multirow{2}{*} \textbf{dparm[4] } & $HIPS\_AMALG$ : Numerical threshold for
coupling between the \\
& interior level and Schur: Default=0.005\\

\multirow{2}{*} \textbf{dparm[5] } & $HIPS\_DROPSCHUR$ : Numerical threshold
for coupling between the \\
& interior level and Schur: Default=0.005\\

 \end{tabular}
\end{center}
\caption{Significations of \textbf{dparams} corresponding to HIPS interface }

\end{table}

### Interfacing with HYPRE
\textbf{HYPRE} ( \textit{High Level Preconditioner}) is a suite of parallel
preconditioner developed at
 Lawrence Livermore National Lab \cite{HYPRE} $\codered$ .

There are two main classes of preconditioners developed in HYPRE: AMG
(Algebraic MultiGrid) and Parasails (Parallel Sparse Approximate
Inverse).

Now, suppose we want to solve $Ax=b$.
At the heart of AMG there is a series of progressively coarser(smaller)
representations of the matrix $A$. Given an approximation
$\hat{x}$ to the solution $x$, consider solving the residual equation $Ae=r$ to
find the error $e$, where $r=b-A\hat{x}$.
A fundamental principle of AMG is that it is an algebraically smooth error. To
reduce the algebraically smooth errors further, they
need to be represented by a smaller defect equation (coarse grid residual
equation) $A_ce_c=r_c$, which is cheaper to solve.
After solving this coarse equation, the solution is then interpolated in fine
grid represented here by matrix $A$. The quality of
AMG depends on the choice of coarsening and interpolating operators.

The \textit{sparse approximate inverse } approximates the inverse of a matrix
$A$ by a sparse matrix $M$. A technical idea to
construct matrix $M$ is to minimize the Frobenuis norm of the residual matrix
$I-MA$. For more details on this preconditioner
technics see \cite{chow} $\codered$.

HYPRE implement three Krylov subspace solvers: GMRES, PCG and BiCGStab.

\paragraph{Installation of HYPRE}
To install HYPRE, first download the HYPRE package at \cite{HYPRE} $\codered$, unpack it
and go to the HYPRE/src source directory and do
\textbf{./configure} to configure Hypre. After this just type \textit{make all}
to create \textbf{libHYPRE.a}.
\paragraph*{Using HYPRE as interface to FreeFem++}
Before calling HYPRE solver inside FreeFem++, you must
compile the file $hypre\_FreeFem.cpp$ to create dynamic library
$hypre\_FreeFem.so$.
To do this, move to the directory $src/solver$ of FreeFem++, edit the file
$makefile.inc$ to specify the following variables:\\
\begin{tabular}{ll}
 %\begin{tabular*}
\textbf{$HYPRE\_DIR$} : & Directory of HYPRE \\
\multirow{2}{*} \textbf{$HYPRE\_INCLUDE$} = &
-I\$($HYPRE\_DIR$)src/hypre/include/ :\\
& Directory for header of HYPRE\\
\textbf{$HYPRE\_LIB$} = & -L\$($HIPS\_DIR$)/src/lib/ -lHYPRE : Hypre Library \\
\textbf{$FREEFEM$} : & FreeFem++ directory \\
\textbf{$FREEFEM\_INCLUDE$} : & FreeFem header for sparse linear solver\\
\textbf{$METIS$} : & METIS directory \\
\textbf{$METIS\_LIB$} : & METIS library \\
\textbf{$MPI$} : & MPI directory \\
\textbf{$MPI\_INCLUDE$} : & MPI headers \\
\end{tabular}


Like with pARMS, the calling of HIPS in FreeFem++ can be done in three manners.
We will present only one example where the user specifies its parameters through
keywords `:::freefem lparams` and `:::freefem dparams`.

\paragraph*{Laplacian 3D solve with HYPRE}
Let us consider again the 3D Laplacian example inside FreeFem++ package where
after discretization we want to solve the linear equation with Hypre. Example
\ref{hypre:laplacian} $\codered$ is the Laplacian3D using Hypre as linear solver.
Example \ref{hypre:laplacian} $\codered$ is the same as \ref{hips:laplacian} $\codered$, so we just
show here the lines where we set some Hypre parameters.

We first load the Hypre solver at line 2. From line 4 to 15 we specifies the
parameters to set to Hypre solver and in line 43
we set parameters to Hypre solver.

 It should be noted that the meaning
of the entries of these vectors is different from those of Hips .
In the case of HYPRE, the meaning of differents entries of vectors
\textbf{iparm} and \textbf{dparm} are given in tables \ref{communipramsHypre} $\codered$ to
\ref{AMGHypre} $\codered$.\\

In Table \ref{hyprelabel} $\codered$ the results of running example \ref{hypre:laplacian} $\codered$
on Cluster Paradent of Grid5000 are reported. We can see in this running
example the efficiency of parallelism, in particular when AMG are use as
preconditioner.
\begin{example}[Laplacian3D.edp]

```freefem
1: load "msh3"
2: load "hipre_FreeFem" //load librairie
3: int nn=10,iii;
4: int[int] iparm(20);
5: real[int] dparm(6);
6: for(iii=0;iii<20;iii++)iparm[iii]=-1;
7: for(iii=0;iii<6;iii++) dparm[iii]=-1;
8: iparm[0]=2; //PCG as krylov method
9: iparm[1]=0; //AMG as preconditionner 2: if ParaSails
10:iparm[7]=7; //Interpolation
11:iparm[9]=6; //AMG Coarsen type
12: iparm[10]=1; //Measure type
13: iparm[16]=2; //Additive schwarz as smoother
13:dparm[0]=1e-13;//Tolerance to convergence
14: dparm[1]=5e-4; //Threshold
15: dparm[2]=5e-4; //truncation factor
.
.
.
43: set(Aa,solver=sparsesolver,dparams=dparm, lparams=iparm);
```

\end{example}



\begin{table}[hbtp]
 \begin{tabular}{|l|l|}
\multirow{2}{*}{iparms[0]} & Solver identification: \\
& 0: BiCGStab, 1: GMRES, 2: PCG. By \textbf{default=1}\\
\multirow{2}{*}{iparms[1]} & Preconditioner identification: \\
& 0: BOOMER AMG, 1: PILUT, 2: Parasails, 3: Schwartz \textbf{Default=0}\\
iparms[2] & Maximum of iteration: \textbf{Default=1000} \\
iparms[3] & Krylov subspace dim: \textbf{Default= 40} \\
iparms[4] & Solver print info level: \textbf{Default=2} \\
iparms[5] & Solver log : \textbf{Default=1} \\
iparms[6] & Solver stopping criteria only for BiCGStab : \textbf{Default=1}
\\
dparms[0] & Tolerance for convergence : \textbf{$Default=1.0e-11$} \\
 \end{tabular}
\caption{Definitions of common entries of \textbf{iparms} and \textbf{dparms}
vectors for every preconditioner in HYPRE}

\end{table}

\begin{table}[hbtp]
 \begin{tabular}{|l|l|}
iparms[7] & AMG interpolation type: \textbf{Default=6} \\
\multirow{2}{*}{iparms[8]} & Specifies the use of GSMG - geometrically \\
& smooth coarsening and interpolation: \textbf{Default=1} \\
iparms[9] & AMG coarsen type: \textbf{Default=6} \\
\multirow{2}{*}{iparms[10]} & Defines whether local or global measures \\
& are used: \textbf{Default=1}\\
iparms[11]& AMG cycle type:\textbf{ Default=1}\\
iparms[12]& AMG Smoother type: \textbf{Default=1}\\
iparms[13]& AMG number of levels for smoothers: \textbf{Default=3}\\
iparms[14]& AMG number of sweeps for smoothers: \textbf{Default=2}\\
iparms[15]& Maximum number of multigrid levels:\textbf{ Default=25}\\
\multirow{6}{*}{iparms[16]}& Defines which variant of the Schwartz method is
used:\\
& 0: hybrid multiplicative Schwartz method (no overlap across processor
boundaries)\\
& 1: hybrid additive Schwartz method (no overlap across processor boundaries)\\
& 2: additive Schwartz method\\
& 3: hybrid multiplicative Schwartz method (with overlap across processor
boundaries)\\
& \textbf{ Default=1}\\
iparms[17]& Size of the system of PDEs: \textbf{ Default=1}\\
iparms[18]& Overlap for the Schwarz method: \textbf{ Default=1}\\
\multirow{4}{*}{iparms[19]} & Type of domain used for the Schwarz method\\
& 0: each point is a domain \\
& 1: each node is a domain (only of interest in ``systems'' AMG)\\
& 2: each domain is generated by agglomeration (default) \\
dparms[1]& AMG strength threshold: \textbf{ Default=0.25}\\
dparms[2]& Truncation factor for the interpolation: \textbf{ Default=1e-2} \\

\multirow{2}{*}{dparms[3]}& Sets a parameter to modify the definition \\
& of strength for diagonal dominant portions of the matrix: \textbf{
Default=0.9} \\
\multirow{2}{*}{dparms[3]} & Defines a smoothing parameter for the additive
Schwartz method \\
& \textbf{ Default=1.} \\
 \end{tabular}
\caption{Definitions of other entries of \textbf{iparms} and \textbf{dparms} if
preconditioner is \textbf{BOOMER AMG}}

\end{table}

\begin{table}[hbtp]
\begin{center}
 \begin{tabular}{|l|l|}
iparms[7]& Row size in Parallel ILUT: \textbf{ Default=1000} \\
iparms[8]& Set maximum number of iterations: \textbf{ Default=30} \\
dparms[1]& Drop tolerance in Parallel ILUT: \textbf{ Default=1e-5} \\
 \end{tabular}
\end{center}
\caption{Definitions of other entries of \textbf{iparms} and \textbf{dparms} if
preconditioner is \textbf{PILUT}}

\end{table}

\begin{table}[hbtp]
 \begin{tabular}{|l|l|}
iparms[7]& Number of levels in Parallel Sparse Approximate inverse: \textbf{
Default=1} \\
\multirow{4}{*}{iparms[8]}& Symmetric parameter for the ParaSails
preconditioner:\\
& 0: nonsymmetric and/or indefinite problem, and nonsymmetric preconditioner\\
& 1: SPD problem, and SPD (factored) preconditioner\\
& 2: nonsymmetric, definite problem, and SPD (factored) preconditioner\\
& \textbf{ Default=0} \\
\multirow{3}{*}{dparms[1]} & Filters parameters:The filter parameter is used to
\\
& drop small nonzeros in the preconditioner, to reduce \\
& the cost of applying the preconditioner: \textbf{ Default=0.1} \\
dparms[2] & Threshold parameter: \textbf{ Default=0.1} \\
 \end{tabular}
\caption{Definitions of other entries of \textbf{iparms} and \textbf{dparms} if
preconditioner is \textbf{ParaSails}}

\end{table}

\begin{table}[hbtp]
 \begin{tabular}{|l|l|}
\multirow{6}{*}{iparms[7]}& Defines which variant of the Schwartz method is
used:\\
& 0: hybrid multiplicative Schwartz method (no overlap across processor
boundaries)\\
& 1: hybrid additive Schwartz method (no overlap across processor boundaries)\\
& 2: additive Schwartz method\\
& 3: hybrid multiplicative Schwartz method (with overlap across processor
boundaries)\\
& \textbf{ Default=1}\\

iparms[8]& Overlap for the Schwartz method: \textbf{ Default=1}\\
\multirow{4}{*}{iparms[9]} & Type of domain used for the Schwartz method\\
& 0: each point is a domain \\
& 1: each node is a domain (only of interest in ``systems'' AMG)\\
& 2: each domain is generated by agglomeration (default) \\
 \end{tabular}
\caption{Definitions of other entries of \textbf{iparms} and \textbf{dparms} if
preconditionner is \textbf{Schwartz}}

\end{table}

\begin{table}
\begin{center}
\begin{tabular}{|c|c|c|}

\multicolumn{1}{|c||}{\textbf{n=$ 4 \times 10^6 $}} &
\multicolumn{1}{|c||}{\textbf{nnz=$13\times10^6$}} &
\multicolumn{1}{|c|}{\textbf{Te=571,29}}\\
\multirow{1}{*}{np} & \multicolumn{2}{|c|}{AMG} \\ \cline{2-3}
& nit & time \\
8 & 6 & 1491.83 \\
16 & 5 & 708.49 \\
32 & 4 & 296.22 \\
64 & 4 & 145.64 \\
\end{tabular}
\end{center}
\caption{Convergence and time for solving linear system from example
\ref{exm:segond} $\codered$ }

\end{table}
### Conclusion
With the different runs presented here, we wanted to illustrate the gain in time
when we increase the number of processors used for the simulations. We saw that in
every case the time for the construction of the finite element matrix is constant. This is normal
because until now this phase is sequential in FreeFem++. In contrast, phases
for solving the linear system are parallel. We saw on several examples
presented here that when we increase the number of processors, in general we
decrease the time used for solving the linear systems. But this not true in every
case. In several case, when we increase the number of processors the time to
convergence also increases. There are two main reasons for this.
First, the increase of processors can lead to the increase of volume of exchanged
data across processors consequently increasing the time for solving the linear
systems.

Furthermore, in decomposition domain type preconditioners, the number of processors
generally corresponds to the number of sub domains. In subdomain methods,
generally when we increase the number of subdomains we decrease convergence
quality of the preconditioner. This can increase the time used for
solving linear equations.

To end this, we should note that good use of the preconditioners interfaced in
FreeFem++ is empiric, because it is difficult to know what is a good
preconditioner for some type of problems. Although, the efficiency of
preconditioners sometimes depends on how its parameters are set. For
this reason we advise the user to pay attention to the meaning of the parameters in the
user guide of the iterative solvers interfaced in FreeFem++.


## Domain decomposition
In the previous section, we saw that the phases to construct a matrix are
sequential. One strategy to construct the matrix in parallel is to divide
geometrically the domain into subdomains. In every subdomain we construct a local
submatrix and after that we assemble every submatrix to form the global matrix.

We can use this technique to solve pde directly in domain $\Omega$. In this case,
in every subdomains you have to define artificial boundary conditions to form
consistent equations in every subdomains. After this, you solve equation in
every subdomains and define a strategy to obtain the global solution.

In terms of parallel programming for FreeFem++, with MPI, this means that the user
must be able to divide processors avaible for computation into subgroups of
processors and also must be able to realize different type of communications in
FreeFem++ script. Here is a wrapper of some MPI functions.
### Communicators and groups
\textbf{Groups}\\
mpiGroup grpe(mpiGroup gp,$KN\_<long>$): Create $MPI\_Group$ from existing group \textbf{gp} by
given vector \\

\textbf{Communicators}\\
Communicators is an abstract MPI object which allows MPI user to communicate
across group of processors. Communicators can be
Intra\-communicators(involves a single group) or Inter\-communicators (involves
two groups). When we not specify type of communicator it will be Intra\-communicators\\



\textbf{mpiComm cc(mpiComm comm, mpiGroup gp)}: Creates a new communicator.
\textit{comm} communicator(handle), \textit{gp} group which is a subset of the
group of \textit{comm} (handle). Return new communicator \\

\textbf{mpiComm cc(mpiGroup gp)}: Same as previous constructor but default
\textit{comm} here is MPI\_COMM\_WORLD.

\textbf{mpiComm cc(mpiComm comm, int color, int key)}: Creates new communicators
based on \textit{colors} and \textit{key}. This constructor is based on
MPI\_Comm\_split routine of MPI.


\textbf{mpiComm cc(MPIrank p,int key)}: Same constructor than the last one.
Here \textit{colors} and \textit{comm} is defined in \textit{MPIrank}. This
constructor is based on MPI\_Comm\_split routine of MPI.
\begin{example}[commsplit.edp]

```freefem
1: int color=mpiRank(comm)\%2;
2: mpiComm ccc(processor(color,comm),0);
3: mpiComm qpp(comm,);
4: mpiComm cp(cc,color,0);
```

\end{example}

\textbf{mpiComm cc(mpiComm comm, int high)}: Creates an intracommunicator from
an intercommunicator. \textit{comm} intercommunicator,
\textit{high} Used to order the groups within \textit{comm} (logical) when
creating the new communicator. This constructor is based on
MPI\_Intercomm\_merge routine of MPI.

\textbf{mpiComm cc(MPIrank p1, MPIrank p2, int tag)}: This constructor creates
an intercommuncator from two intracommunicators. \textit{p1} defined local
(intra)communicator and rank in local\_comm of leader (often 0) while
\textit{p2} defined remote communicator and rank in peer\_comm of remote leader
(often 0). \textit{tag} Message tag to use in constructing intercommunicator.
This constructor is based on MPI\_Intercomm\_create.

\begin{example}[merge.edp]

```freefem
1: mpiComm comm,cc;
2: int color=mpiRank(comm)\%2;
3: int rk=mpiRank(comm);
4: int size=mpiSize(comm);
4: cout << "Color values " << color << endl;
5: mpiComm ccc(processor((rk<size/2),comm),rk);
6: mpiComm cp(cc,color,0);
7: int rleader;
8: if (rk == 0) { rleader = size/2; }
9: else if (rk == size/2) { rleader = 0;}
10: else { rleader = 3; }
11: mpiComm qqp(processor(0,ccc),processor(rleader,comm),12345);
12:int aaa=mpiSize(ccc);
13:cout << "number of processor" << aaa << endl;
```

\end{example}

### Process
In FreeFem++ we wrap MPI process by function call \textbf{processor} which
create internal FreeFem++ object call \textbf{MPIrank}. This mean that do not
use \textbf{MPIrank} in FreeFem++ script.

\textbf{processor(int rk):} Keep process rank inside object \textbf{MPIrank}.
Rank is inside MPI\_COMM\_WORLD.

\textbf{processor(int rk, mpiComm cc) and processor(mpiComm cc,int rk) } process
rank inside communicator cc.

\textbf{processor(int rk, mpiComm cc) and processor(mpiComm cc,int rk) } process
rank inside communicator cc.

\textbf{processorblock(int rk) }: This function is exactlly the same than
\textbf{processor(int rk)} but is use in case of blocking communication.

\textbf{processorblock(int rk, mpiComm cc) }: This function is exactlly the same
than \textbf{processor(int rk,mpiComm cc)} but use a synchronization point.



### Points to Points communicators
In FreeFem++ you can call MPI points to points communications functions.



\textbf{Send(processor(int rk,mpiComm cc),Data D)} : Blocking send of
\textit{Data D} to processor of \textit{rank rk} inside communicator
\textit{cc}. Note that \textit{Data D} can be: \textit{int, real,complex ,
int[int], real[int],complex[int], Mesh, Mesh3, Matrix}.


\textbf{Recv(processor(int rk,mpiComm cc),Data D)}: Receive \textit{Data D} from
process of rank \textit{rk} in communicator \textit{cc}. Note that \textit{Data
D} can be: \textit{int, real,complex , int[int], real[int],complex[int], Mesh,
Mesh3, Matrix} and should be the same type than corresponding send.


\textbf{Isend(processor(int rk,mpiComm cc),Data D)} : Non blocking send of
\textit{Data D} to processor of \textit{rank rk} inside communicator
\textit{cc}. Note that \textit{Data D} can be: \textit{int, real,complex ,
int[int], real[int],complex[int], Mesh, Mesh3, Matrix}.

\textbf{Recv(processor(int rk,mpiComm cc),Data D)}: Receive corresponding to
send.

### Global operations

In FreeFem++ you can call MPI global communication functions.

\textbf{broadcast(processor(int rk,mpiComm cc),Data D)}: Process \textit{rk}
Broadcast \textit{Data D} to all process inside \textit{communicator cc}. Note
that \textit{Data D} can be: \textit{int, real,complex , int[int],
real[int],complex[int], Mesh, Mesh3, Matrix}.\\


\textbf{broadcast(processor(int rk),Data D)}: Process \textit{rk} Broadcast
\textit{Data D} to all process inside

MPI\_COMM\_WORLD. Note that \textit{Data D} can be: \textit{int, real,complex ,
int[int], real[int],complex[int], Mesh, Mesh3, Matrix}.\\


\textbf{mpiAlltoall(Data a,Data b)}: Sends \textit{data a} from all to all
processes. Receive buffer is \textit{Data b}. This is done inside communicator
MPI\_COMM\_WORLD.\\


\textbf{mpiAlltoall(Data a,Data b, mpiComm cc)}: Sends \textit{data a} from all
to all processes. Receive buffer is \textit{Data b}. This is done inside
communicator cc.\\


\textbf{mpiGather(Data a,Data b,processor(mpiComm,int rk)} : Gathers together
values \textit{Data a} from a group of processes. Process of rank \textit{rk}
get data on communicator \textit{rk}. This function is like MPI\_Gather\\


\textbf{mpiAllgather(Data a,Data b)} : Gathers \textit{Data a} from all
processes and distribute it to all in \textit{Data b}. This is done inside
communicator MPI\_COMM\_WORLD. This function is like MPI\_Allgather\\

\textbf{mpiAllgather(Data a,Data b, mpiComm cc)} : Gathers \textit{Data a} from
all processes and distribute it to all in \textit{Data b}. This is done inside
\textbf{communicator cc}. This function is like MPI\_Allgather\\



\textbf{mpiScatter(Data a,Data b,processor(int rk, mpiComm cc))} : Sends
\textbf{Data a} from one process whith rank \textbf{rk} to all other processes
in group represented by communicator \textit{mpiComm cc}.\\

\textbf{mpiReduce(Data a,Data b,processor(int rk, mpiComm cc),MPI\_Op op), }
Reduces values \textit{Data a} on all processes
to a single value \textit{Data b} on process of rank \textit{rk} and
communicator \textit{cc}.
Operation use in reduce is: \textit{MPI\_Op op} which can be: \textit{mpiMAX},
\textit{mpiMIN}, \textit{mpiSUM},
 \textit{mpiPROD}, \textit{mpiLAND}, \textit{mpiLOR}, \textit{mpiLXOR},
\textit{mpiBAND},
 \textit{mpiBXOR}, \textit{mpiMAXLOC}, \textit{mpiMINLOC}.

Note that, for all global operations, only int[int] and real[int] are data type
take in account in FreeFem++.

% exemple obsolete F. Hecht revmove 10/08/2016
% ---------------
%The following example present in details of Schwartz domain decomposition
%algorithm for solving Laplacian2d problem. In this example we use two level
%of parallelism to solve simple Laplacian2d in square domain. We have few number
%of subdomain and in every subdomain we use parallel sparse solver to solve local problem.
%
%
%\begin{example}[schwarz.edp]
%
```freefem
%//
%1:load "hypre_FreeFem"; //Load Hypre solver
%2:func bool AddLayers(mesh & Th,real[int] &ssd,int n,real[int] &unssd)
%{
% // build a continuous function uussd (P1) :
% // ssd in the caracteristics function on the input sub domain.
% // such that :
% // unssd = 1 when ssd =1;
% // add n layer of element (size of the overlap)
% // and unssd = 0 ouside of this layer ...
% // ---------------------------------
% fespace Vh(Th,P1);
% fespace Ph(Th,P0);
% Ph s;
% assert(ssd.n==Ph.ndof);
% assert(unssd.n==Vh.ndof);
% unssd=0;
% s[]= ssd;
% // plot(s,wait=1,fill=1);
% Vh u;
% varf vM(u,v)=int2d(Th,qforder=1)(u*v/area);
% matrix M=vM(Ph,Vh);
%
% for(int i=0;i<n;++i)
% {
% u[]= M*s[];
% // plot(u,wait=1);
% u = u>.1;
% // plot(u,wait=1);
% unssd+= u[];
% s[]= M'*u[];//';
% s = s >0.1;
% }
% unssd /= (n);
% u[]=unssd;
% ssd=s[];
% return true;
%}
%3: mpiComm myComm; // Create communicator with value MPI\_COMM\_WORLD
%
%4: int membershipKey,rank,size; //Variables for manage communicators
%5: rank=mpiRank(myComm); size=mpiSize(myComm); //Rank of process and size of communicator
%6: bool withmetis=1, RAS=0; //Use or not metis for partitioning Mesh
%7: int sizeoverlaps=5; // size off overlap
%8: int withplot=1;
%9: mesh Th=square(100,100);
%10: int[int] chlab=[1,1 ,2,1 ,3,1 ,4,1 ];
%11: Th=change(Th,refe=chlab);
%12: int nn=2,mm=2, npart= nn*mm;
%13: membershipKey = mpiRank(myComm)\%npart; // Coloring for partitioning process group
%14: mpiComm cc(processor(membershipKey,myComm),rank); //Create MPI communicator according previous coloring
%15: fespace Ph(Th,P0),fespace Vh(Th,P1);
%16: Ph part;
%17: Vh sun=0,unssd=0;
%18: real[int] vsum=sun[],reducesum=sun[]; //Data use for control partitioning.
%19: Ph xx=x,yy=y,nupp;
%20: part = int(xx*nn)*mm + int(yy*mm);
%21: if(withmetis)
% {
% load "metis";
% int[int] nupart(Th.nt);
% metisdual(nupart,Th,npart);
% for(int i=0;i<nupart.n;++i)
% part[][i]=nupart[i];
% }
%22: if(withplot>1)
%21: plot(part,fill=1,cmm="dual",wait=1);
%22: mesh[int] aTh(npart);
%23: mesh Thi=Th;
%24: fespace Vhi(Thi,P1);
%25: Vhi[int] au(npart),pun(npart);
%26: matrix[int] Rih(npart), Dih(npart), aA(npart);
%27: Vhi[int] auntgv(npart), rhsi(npart);
%28: i=membershipKey;
% Ph suppi= abs(part-i)<0.1;
% AddLayers(Th,suppi[],sizeoverlaps,unssd[]);
% Thi=aTh[i]=trunc(Th,suppi>0,label=10,split=1);
% Rih[i]=interpolate(Vhi,Vh,inside=1); // Vh -> Vhi
% if(RAS)
% {
% suppi= abs(part-i)<0.1;
% varf vSuppi(u,v)=int2d(Th,qforder=1)(suppi*v/area);
% unssd[]= vSuppi(0,Vh);
% unssd = unssd>0.;
% if(withplot>19)
% plot(unssd,wait=1);
% }
% pun[i][]=Rih[i]*unssd[];//this is global operation
% sun[] += Rih[i]'*pun[i][];// also global operation like broadcast';
% vsum=sun[];
% if(withplot>9)
% plot(part,aTh[i],fill=1,wait=1);
% // Add mpireduce for sum all sun and pun local contribution.
%29: mpiReduce(vsum, reducesum,processor(0,myComm),mpiSUM); //MPI global operation MPi\_Reduce on global communicator
%30: broadcast(processor(0,myComm),reducesum); //Broadcast sum on process 0 to all process
%31: sun[]=reducesum;
%32: plot(sun,wait=1);
%33: i=membershipKey
%34: Thi=aTh[i];
%35: pun[i]= pun[i]/sun;
%36: if(withplot>8) plot(pun[i],wait=1);
%37: macro Grad(u) [dx(u),dy(u)]//EOM
%38: sun=0;
%39: i=membershipKey
% Thi=aTh[i];
% varf va(u,v) =
% int2d(Thi)(Grad(u)'*Grad(v))//')
% +on(1,u=1) + int2d(Th)(v)
% +on(10,u=0) ;
%40: aA[i]=va(Vhi,Vhi);
%41: set(aA[i],solver=sparsesolver,mpicomm=cc); //Set parameters for Solver Hypre. mpicomm=cc means you not solve on global process but in group on of process define by cc
%42: rhsi[i][]= va(0,Vhi);
%43: Dih[i]=pun[i][];
%44: real[int] un(Vhi.ndof);
%45: un=1.;
%46: real[int] ui=Dih[i]*un;
%47: sun[] += Rih[i]'*ui;//';
%48: varf vaun(u,v) = on(10,u=1);
%49: auntgv[i][]=vaun(0,Vhi); // store arry of tgv on Gamma intern.
%56: reducesum=0; vsum=sun;
%57: mpiReduce(vsum, reducesum,processor(0,myComm),mpiSUM); //MPI global operation MPi\_Reduce on global communicator
%58: broadcast(processor(0,myComm),reducesum); //Broadcast sum on process 0 to all other process
%59: sun[]=reducesum;
%60: if(withplot>5)
%61: plot(sun,fill=1,wait=1);
%62: cout << sun[].max << " " << sun[].min<< endl;
%63: assert( 1.-1e-9 <= sun[].min && 1.+1e-9 >= sun[].max);
%64: int nitermax=1000;
%{
% Vh un=0;
% for(int iter=0;iter<nitermax;++iter)
% {
% real err=0,rerr=0;
% Vh un1=0;
% i=membershipKey;
% Thi=aTh[i];
% real[int] ui=Rih[i]*un[];//';
%	 real[int] bi = ui .* auntgv[i][];
% bi = auntgv[i][] ? bi : rhsi[i][];
% ui=au[i][];
% ui= aA[i] ^-1 * bi; //Solve local linear system on group of process represented by color membershipKey
% bi = ui-au[i][];
% err += bi'*bi;//';
%	 au[i][]= ui;
% bi = Dih[i]*ui; //Prolongation of current solution to obtain right hand
% un1[] += Rih[i]'*bi;// ';
%}
%65: reducesum=0; vsum=un1[];
%66: mpiReduce(vsum, reducesum,processor(0,myComm),mpiSUM); //MPI global operation MPi\_Reduce on global communicator
%67: broadcast(processor(0,myComm),reducesum); //Broadcast sum on process 0 to all other process
%68: un1[]=reducesum;
%69: real residrela=0;
%70: mpiReduce(err,residrela ,processor(0,myComm),mpiSUM);
%71: broadcast(processor(0,myComm),residrela);
%72: err=residrela; err= sqrt(err);
%73: if(rank==0)	cout << iter << " Err = " << err << endl;
%74: if(err<1e-5) break;
%75: un[]=un1[];
%76: if(withplot>2)
%77: plot(au,dim=3,wait=0,cmm=" iter "+iter,fill=1 );
%78: }
%79: plot(un,wait=1,dim=3);
%80: }
%```

%\end{example}
%
%
%
% ------





\input{petschpddm}
