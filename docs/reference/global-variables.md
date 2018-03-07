## CG

Conjugate gradient solver.

> Usable in [`:::freefem problem`](types/#problem) and [`:::freefem solve`](types/#solve) definition
```freefem
problem Laplacian (U, V, solver=CG) = ...
```

> Or in [`:::freefem matrix`](types/#matrix) construction
```freefem
matrix A = vLaplacian(Uh, Uh, solver=CG);
```

> Or in [`:::freefem set` function](functions/#set)
```freefem
set(A, solver=CG);
```

## Cholesky

Cholesky solver.

## Crout

Crout solver.

## GMRES

GMRES solver (Generalized minimal residual method).

## include
Include a module.
```freefem
include "iovtk"
```

## load
Load a script.
```freefem
load "getARGV.idp"
```

## LU

LU solver.

## N
Normal.
```freefem
func Nx = N.x;
func Ny = N.y;
func Nz = N.z;
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

## sparsesolver

Sparse matrix solver.

## verbosity
Verbosity level.
```freefem
int Verbosity = verbosity;
verbosity = 0;
```
0 = nothing, 1 = few informations, 10 = a lot of informations, ...

This is an integer value.

## version
FreeFem++ version.
```freefem
cout << version << endl;
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
real CurrentY = y;
```
This is a real value.
