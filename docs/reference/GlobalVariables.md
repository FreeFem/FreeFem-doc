## area

Area of the current triangle.

```freefem
fespace Vh0(Th, P0);
Vh0 A = area;
```

## ARGV

Array that contains all the command line arguments.

```freefem
for (int i = 0; i < ARGV.n; i++)
	cout << ARGV[i] << endl;
```

See [Command line arguments example](/examples/#command-line-arguments) for a complete example.

## BoundaryEdge

Return 1 if the current edge is on a boundary, 0 otherwise.

```freefem
real B = int2d(Th)(BoundaryEdge);
```

## CG

Conjugate gradient solver.

> Usable in [`:::freefem problem`](../Types/#problem) and [`:::freefem solve`](../Types/#solve) definition
```freefem
problem Laplacian (U, V, solver=CG) = ...
```

> Or in [`:::freefem matrix`](../Types/#matrix) construction
```freefem
matrix A = vLaplacian(Uh, Uh, solver=CG);
```

> Or in [`:::freefem set` function](../Functions/#set)
```freefem
set(A, solver=CG);
```

## Cholesky

Cholesky solver.

## Crout

Crout solver.

## edgeOrientation
Sign of $i-j$ if the current edge is $[q_i, q_j]$.

```freefem
real S = int1d(Th, 1)(edgeOrientation);
```

## false

False boolean value.

```freefem
bool b = false;
```

## GMRES

GMRES solver (Generalized minimal residual method).

## hTriangle

Size of the current triangle.

```freefem
fespace Vh(Th, P0);
Vh h = hTriangle;
```

## include
Include an [external library](../ExternalLibraries).
```freefem
include "iovtk"
```

## InternalEdge

Return 0 if the current edge is on a boundary, 1 otherwise.

```freefem
real I = int2d(Th)(InternalEdge);
```

## label
Label number of a boundary if the current point is on a boundary, 0 otherwise.

```freefem
int L = Th(xB, yB).label;
```

## lenEdge
Length of the current edge.

For an edge $[q_i, g_j]$, return $|q_i-q_j|$.

```freefem
real L = int1d(Th, 1)(lenEdge);
```

## load
Load a script.
```freefem
load "Element_P3"
```

## LU

LU solver.

## N
Outward unit normal at the current point if it is on a curve defined by a border. `:::freefem N.x, N.y, N.z` are respectively the $x$, $y$ and $z$ components of the normal.

```freefem
func Nx = N.x;
func Ny = N.y;
func Nz = N.z;
```

## nTonEdge

Number of adjacent triangles of the current edge.

```freefem
real nTE = int2d(Th)(nTonEdge);
```

## nuEdge

Index of the current edge in the triangle.

```freefem
real nE = int2d(Th)(nuEdge);
```

## nuTriangle

Index of the current triangle.

```freefem
fespace Vh(Th, P0);
Vh n = nuTriangle;
```

## P
Current point.
```freefem
real cx = P.x;
real cy = P.y;
real cz = P.z;
```

## pi
Pi = 3.14159.
```freefem
real Pi = pi;
```
This is a real value.

## region
Region number of the current point. If the point is outside, then `:::freefem region == notaregion` where `:::freefem notaregion` is a __`FreeFem++`__ integer constant.

```freefem
int R = Th(xR, yR).region;
```

## sparsesolver

Sparse matrix solver.

## true

True boolean value.

```freefem
bool b = true;
```

## verbosity
Verbosity level.
```freefem
int Verbosity = verbosity;
verbosity = 0;
```
0 = nothing, 1 = little information, 10 = a lot of information, ...

This is an integer value.

## version
FreeFem++ version.
```freefem
cout << version << endl;
```

## volume

Volume of the current tetrahedra.

```freefem
fespace Vh0(Th, P0);
Vh0 V = volume;
```

## x
The $x$ coordinate at the current point.
```freefem
real CurrentX = x;
```
This is a real value.

## y
The $y$ coordinate at the current point.
```freefem
real CurrentY = y;
```
This is a real value.

## z
The $z$ coordinate at the current point.
```freefem
real CurrentZ = z;
```
This is a real value.
