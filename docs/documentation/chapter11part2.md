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
$src/solver$ type `:::freefem make hips` to create $hips\_FreeFem.so$. Like with pARMS, the calling of HIPS in FreeFem++ can be done in three different manners. We will present only one example where the user specifies the parameters through keywords `:::freefem lparams` and `:::freefem dparams`.

__Laplacian 3D solve with HIPS__

Let us consider the 3D Laplacian example inside FreeFem++ package where after discretization we want to solve the linear equation with Hips. Example \ref{hips:laplacian} $\codered$ is Laplacian3D using Hips as linear solver. We first load Hips solver at line 2. From line 4 to 15 we specify the parameters for the Hips solver and in line 46 of example \ref{hips:laplacian} $\codered$ we set these parameters in the linear solver.

In Table 1.10 \ref{hipslabel} $\codered$ results of running example 1.1 5\ref{hips:laplacian} $\codered$ on Cluster Paradent of Grid5000 are reported. We can see in this running example the efficiency of parallelism.

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

<table>
	<thead>
		<tr>
			<th colspan="3">Table 11.10 : Legend of table 11.10 are give in table 11.9 \ref{legtableparm} $\codered$</th>
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

<table>
	<thead>
		<tr>
			<th colspan="2">Table 11.11 : Significations of lparams corresponding to HIPS interface</th>
		</tr>
	</thead>
	<tbody>
	    <tr>
	        <td style="font-weight: bold">Entries of iparm</td>
	        <td style="font-weight: bold">Significations of each entries</td>
	    </tr>
	    <tr>
	        <td>iparm[0]</td>
	        <td>Strategy use for solving (Iterative=0 or Hybrid=1 or Direct=2). Defaults values are : Iterative</td>
	    </tr>
	    <tr>
	        <td>iparm[1]</td>
	        <td>Krylov methods. If iparm[0]=0, give type of Krylov methods: 0 for GMRES, 1 for PCG</td>
	    </tr>
	    <tr>
	        <td>iparm[2]</td>
	        <td>Maximum of iterations in outer iteration: default value 1000</td>
	    </tr>
	    <tr>
	        <td>iparm[3]</td>
	        <td>Krylov subspace dimension in outer iteration: default value 40</td>
	    </tr>
	    <tr>
	        <td>iparm[4]</td>
	        <td>Symmetric(=0 for symmetric) and 1 for unsymmetricmatrix: default value 1 (unsymmetric matrix)</td>
	    </tr>
	    <tr>
	        <td>iparm[5]</td>
	        <td>Pattern of matrix are symmetric or not: default value 0</td>
	    </tr>
	    <tr>
	        <td>iparm[6]</td>
	        <td>Partition type of input matrix: default value 0</td>
	    </tr>
	    <tr>
	        <td>iparm[7]</td>
	        <td>Number of level that use the HIPS locally consistentfill-in: Default value 2</td>
	    </tr>
	    <tr>
	        <td>iparm[8]</td>
	        <td>Numbering in indices array will start at 0 or 1: Default value 0</td>
	    </tr>

	    <tr>
	        <td>iparm[9]</td>
	        <td>Scale matrix. Default value 1</td>
	    </tr>
	    <tr>
	        <td>iparm[10]</td>
	        <td>Reordering use inside subdomains for reducingfill-in: Only use for iterative. Default value 1</td>
	    </tr>
	    <tr>
	        <td>iparm[11]</td>
	        <td>Number of unknowns per node in the matrix non-zeropattern graph: Default value 1</td>
	    </tr>
	    <tr>
	        <td>iparm[12]</td>
	        <td>This value is used to set the number of time the normalization is applied to the matrix: Default 2.</td>
	    </tr>
	    <tr>
	        <td>iparm[13]</td>
	        <td>Level of informations printed during solving: Default 5.</td>
	    </tr>
	    <tr>
	        <td>iparm[14]</td>
	        <td>HIPS_DOMSIZE Subdomain size</td>
	    </tr>
	</tbody>
</table>

<table>
	<thead>
		<tr>
			<th colspan="2">Table 11.12 : Significations of dparams corresponding to HIPS interface</th>
		</tr>
	</thead>
	<tbody>
	    <tr>
	        <td>dparm[0]</td>
	        <td>$HIPS\_PREC$: Relative residual norm: Default=1e-9</td>
	    </tr>
	    <tr>
	        <td>dparm[1]</td>
	        <td>$HIPS\_DROPTOL0$: Numerical threshold in ILUT for interior domain (important : set 0.0 in HYBRID: Default=0.005)</td>
	    </tr>
	    <tr>
	        <td>dparm[2]</td>
	        <td>$HIPS\_DROPTOL1$ : Numerical threshold in ILUT for Schur preconditioner: Default=0.005</td>
	    </tr>
	    <tr>
	        <td></td>
	        <td></td>
	    </tr>
	    <tr>
	        <td>dparm[3]</td>
	        <td>$HIPS\_DROPTOLE$ : Numerical threshold for coupling between the interior level and Schur: Default 0.005</td>
	    </tr>
	    <tr>
	        <td>dparm[4]</td>
	        <td>$HIPS\_AMALG$ : Numerical threshold for coupling between the interior level and Schur: Default=0.005</td>
	    </tr>
	    <tr>
	        <td>dparm[5]</td>
	        <td>$HIPS\_DROPSCHUR$ : Numerical threshold for coupling between the interior level and Schur: Default=0.005</td>
	    </tr>
	</tbody>
</table>

### Interfacing with HYPRE

__HYPRE__ (High Level Preconditioner) is a suite of parallel preconditioner developed at Lawrence Livermore National Lab \cite{HYPRE} $\codered$ .

There are two main classes of preconditioners developed in HYPRE: AMG (Algebraic MultiGrid) and Parasails (Parallel Sparse Approximate Inverse).

Now, suppose we want to solve $Ax=b$.
At the heart of AMG there is a series of progressively coarser(smaller) representations of the matrix $A$. Given an approximation $\hat{x}$ to the solution $x$, consider solving the residual equation $Ae=r$ to find the error $e$, where $r=b-A\hat{x}$. A fundamental principle of AMG is that it is an algebraically smooth error. To reduce the algebraically smooth errors further, they need to be represented by a smaller defect equation (coarse grid residual equation) $A_ce_c=r_c$, which is cheaper to solve. After solving this coarse equation, the solution is then interpolated in fine grid represented here by matrix $A$. The quality of AMG depends on the choice of coarsening and interpolating operators.

The _sparse approximate inverse_ approximates the inverse of a matrix $A$ by a sparse matrix $M$. A technical idea to construct matrix $M$ is to minimize the Frobenuis norm of the residual matrix $I-MA$. For more details on this preconditioner technics see \cite{chow} $\codered$.

HYPRE implement three Krylov subspace solvers: GMRES, PCG and BiCGStab.

__Installation of HYPRE__

To install HYPRE, first download the HYPRE package at \cite{HYPRE} $\codered$, unpack it and go to the HYPRE/src source directory and do `:::bash ./configure` to configure Hypre. After this just type `:::bash make all` to create __libHYPRE.a__.

__Using HYPRE as interface to FreeFem++__

Before calling HYPRE solver inside FreeFem++, you must compile the file $hypre\_FreeFem.cpp$ to create dynamic library $hypre\_FreeFem.so$.
To do this, move to the directory $src/solver$ of FreeFem++, edit the file $makefile.inc$ to specify the following variables:

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
