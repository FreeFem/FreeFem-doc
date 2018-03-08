## Standard types

### int
Integer value.
```freefem
int i = 0;
```

### bool
Boolean value.
```freefem
bool b = true;
```

### real
Real value (C double precision).
```freefem
real r = 0.;
```

### complex
Complex value (two C double precision).
```freefem
complex c = 0. + 1i;
```
The imaginary number $i$ is defined as `1i`

### string
String value.
```freefem
string s = "this is a string";
```


## Mesh design

### border
Border type.
```freefem
border b(t=0., 1.){x=cos(2.*pi*t); y=sin(2.*pi*t);}
```
Define the 2D geometrical border in parametric coordinates.

!!!note "Label"
	A label can be defined with the border:
	```freefem
	border b(t=0., 1.){x=cos(2.*pi*t); y=sin(2.*pi*t); label=1;}
	```

!!!note "Inner variable"
	An inner variable can be defined inside a border:
	```freefem
	border b(t=0., 1.){real tt=2.*pi*t; x=cos(tt); y=sin(tt);}
	```

!!!note "From vector"
	A border can be defined from two vectors using `P.x` and `P.y`:
	```freefem
	border b(t=0, vectorX.n-1){x=vectorX[t]; P.x=vectorY[t];}
	```

### mesh
2D Mesh type
```freefem
mesh Th;
```

### mesh3
3D mesh type
```freefem
mesh3 Th;
```


## Finite element space design

### fespace
Finite element space type.
```freefem
fespace Uh(Th, P1);
fespace UPh(Th, [P2, P2, P1]);
```
A finite element space is based on a mesh (`Th`) with an element definition, scalar (`:::freefem P1`) or vector (`:::freefem [P2, P2, P1]`).


**Available finite element space:**

Generic:

 - `:::freefem P0 / P03d`
 - `:::freefem P0Edge`
 - `:::freefem P1 / P13d`
 - `:::freefem P1dc`
 - `:::freefem P1b / P1b3d`
 - `:::freefem P1bl / P1bl3d`
 - `:::freefem P1nc`
 - `:::freefem P2 / P23d`
 - `:::freefem P2b`
 - `:::freefem P2dc`
 - `:::freefem RT0 / RT03d`
 - `:::freefem RT0Ortho`
 - `:::freefem Edge03d`

Using _Element_P3_:

 - `:::freefem P3`

Using _Element_P3dc_:

 - `:::freefem P3dc`

Using _Element_P4_:

 - `:::freefem P4`

Using _Element_P4dc_:

 - `:::freefem P4dc`

Using _Element_PkEdge_:

 - `:::freefem P1Edge`
 - `:::freefem P2Edge`
 - `:::freefem P3Edge`
 - `:::freefem P4Edge`
 - `:::freefem P5Edge`

Using _Morlay_:

 - `:::freefem P2Morley`

Using _HCT_:

 - `:::freefem HCT`

Using _BernardiRaugel_:

 - `:::freefem P2BR`

Using _Element_Mixte_:

 - `:::freefem RT1`
 - `:::freefem RT1Ortho`
 - `:::freefem BDM1`
 - `:::freefem BDM1Ortho`

Using _Element_Mixte3d_:

 - `:::freefem Edge13D`
 - `:::freefem Edge23D`

Using _Element_QF_:

 - `:::freefem FEQF`

A finite element variable is defined as follow:
```freefem
fespace Uh(Th, P1);
Uh u;

fespace UPh(Th, [P2, P2, P1]);
UPh [Ux, Uy, p];
```

## Macro design

### macro
Macro type.
```freefem
macro grad(u) [dx(u), dy(u)] //
```
Macro ends with `//`.

!!!note "Macro concatenation"
	You can use the C concatenation operator ## inside a macro using #.

	If `Ux` and `Uy` are defined as finite element variable, you can define:
	```freefem
	macro Grad(U) [grad(U#x), grad(U#y)] //
	```


## Functions design

### func
Function type.

Function without parameters ($x$, $y$ and $z$ are implicitly considered):
```freefem
func f = x^2 + y^2;
```

Function with parameters:
```freefem
func real f (real var){
	return x^2 + y^2 + var^2;
}
```


## Problem design

### problem
Problem type.
```freefem
problem Laplacian (u, uh) = ...
```
FreeFem++ needs the variational form in the problem definition.

In order to solve the problem, just call:
```freefem
Laplacian;
```

!!!note "Solver"
	A solver can be specified in the problem definition:
	```freefem
	problem Laplacian(u, uh, solver=CG) = ...
	```

	The default solver is `:::freefem sparsesolver` or `:::freefem LU` if any direct sparse solver is available.

	Solvers are listed in the [Global variables](global-variables) section.

!!!note "Stop test"
	A criterion $\varepsilon$ can be defined for iterative methods, like CG for example:
	```freefem
	problem Laplacian(u, uh, solver=CG, eps=1.e-6) = ...
	```

	If $\varepsilon>0$, the stop test is:
	$$
	||Ax-b|| < \varepsilon
	$$
	Else, the stop test is:
	$$
	||Ax-b|| < \frac{|\varepsilon|}{||Ax_0-b||}
	$$

!!!note "Reconstruction"
	The keyword `:::freefem init` controls the reconstruction of the internal problem matrix.

	If `:::freefem init` is set to `:::freefem false` or `:::freefem 0`, the matrix is reconstructed et each problem calls (or after a mesh modification), else the previously constructed matrix is used.

	```freefem
	problem Laplacian(u, uh, init=1) = ...
	```

!!!note "Preconditioning"
	A preconditioner can be specified in the problem definition:

	```freefem
	problem Laplacian(u, uh, precon=P) = ...
	```

	The preconditioning function must have a prototype like:

	```freefem
	func real[int] P(real[int] &xx);
	```


!!!note "_Très grande valeur_"
	The "_Très grand valeur_" `:::freefem tgv` (or Terrible giant value) used to implement the Dirichlet conditions can be modified in the problem definition:

	```freefem
	problem Laplacian(u, uh, tgv=1e30) = ...
	```

	Refere to $\codered$ for a description of the Dirichlet condition implementation.

!!!note "Pivot tolerance"
	The tolerance of the pivot in `:::freefem UMFPACK`, `:::freefem LU`, `:::freefem Crout`, `:::freefem Cholesky` factorization can be modified in the problem definition:

	```freefem
	problem Laplacian(u, uh, solver=LU, tolpivot=1e-20) = ...
	```

!!!note "`:::freefem UMFPACK`"
	Two specific parameters for the `:::freefem UMFPACK` can be modifed:

	 * Tolerance of the pivot sym
	 * strategy

	```freefem
 	problem Laplacian(u, uh, solver=LU, tolpivotsym=1e-1, strategy=0) = ...
 	```

	Refer to the [UMFPACK website](http://faculty.cse.tamu.edu/davis/research.html) for more informations.

Usage of `:::freefem problem` is detailled in the [tutorial](../tutorial).

### solve
Solve type.

Identical to [problem](#problem) but automatically solved.

Usage of `:::freefem solve` is detailled in the [tutorial](../tutorial).

### varf
Variational form type.

```freefem
varf vLaplacian (u, uh) = ...
```

Directly define a variationnal form.

This is the other way to define a problem in order to directly manage matrix and right hang side.

Usage of `:::freefem varf` is detailled in the [tutorial](../tutorial).

## Array

Array can be defined using types: `:::freefem int, bool, real, complex, string, ...`

### Array index
Array index can be int or string:
```freefem
real[int] Ai = [1, 1, 0, 0];
real[string] As = [1, 1, 0, 0];
```

### Array size
The size of an array is obtained using the keyword `n`:
```freefem
int ArraySize = Ai.n;
```

### Double array
A double array (matrix) can be defined using two indexes:

```freefem
real[int, int] Aii = [[1, 1], [0, 0]];
```

The two sizes are obtained using the keywords `n` and `m`:

```freefem
int ArraySize1 = Aii.n;
int ArraySize2 = Aii.m;
```

The minimum and maximum values of an array (simple or double) can be obtained using:
```freefem
real ArrayMin = Aii.min;
real ArrayMax = Aii.max;
```

!!!tip
	An array can be obtained from a finite element variable using:

	```freefem
	real[int] aU = U[];
	```

	where `:::freefem U` is a finite element variable.


## Matrix
Matrices can be defined like vectors:

```freefem
matrix A = [[1, 2, 3],
			[4, 5, 6],
			[7, 8, 9]];
```

or using a variational form type:

```freefem
matrix Laplacian = vLaplacian(Uh, Uh);
```

Matrices are designed using templates, so they can be real or complex:

```freefem
matrix<real> A = ...
matrix<complex> Ai = ...
```

The size of a matrix is obtain using:

```freefem
int NRows = A.n;
int NColumns = A.m;
```

The diagonal of the matrix is obtained using:

```freefem
real[int] Aii = A.diag;
```

!!!note "Solver"
	See [`:::freefem problem`](#problem).

	The default solver is [`:::freefem GMRES`](../global-variables/#GMRES).

	```freefem
	matrix A = vLaplacian(Uh, Uh, solver=sparsesolver);
	```
	or
	```freefem
	set(A , solver=sparsesolver);
	```

!!!note "Factorize"
	If `:::freefem true`, the factorization is done for `:::freefem LU`, `:::freefem Cholesky` or `:::freefem Crout`.

	```freefem
	matrix A = vLaplacian(Uh, Uh, solver=LU, factorize=1);
	```
	or
	```freefem
	set(A , solver=LU, factorize=1);
	```

!!!note "Stop test"
	See [`:::freefem problem`](#problem).

!!!note "_Très grande valeur_"
	See [`:::freefem problem`](#problem).

!!!note "Preconditioning"
	See [`:::freefem problem`](#problem).

!!!note "Pivot tolerance"
	See [`:::freefem problem`](#problem).

!!!note "`:::freefem UMFPACK`"
	See [`:::freefem problem`](#problem).

!!!note "datafilename"
	$\codered$

!!!note "lparams"
	$\codered$

!!!note "dparams"
	$\codered$

!!!note "sparams"
	$\codered$

!!!tip
	To modify the solver, the stop test,... after the matrix construction, use the [`:::freefem set` keyword](functions/#set).
