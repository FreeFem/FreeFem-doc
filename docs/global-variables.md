## CG

Conjugate gradient solver.

> Usable in problem and solve definition
```ffpp
problem Laplacian (U, V, solver=CG) = ...
```
> Or in matrix construction
```ffpp
matrix A = vLaplacian(Uh, Uh, solver=CG);
```

## Cholesky

Cholesky solver.

## Crout

Crout solver.

## GMRES

GMRES solver (Generalized minimal residual method).

## include
Include a module.
```ffpp
include "iovtk"
```

## load
Load a script.
```ffpp
load "getARGV.idp"
```

## LU

LU solver.

## N
Normal.
```ffpp
func Nx = N.x;
func Ny = N.y;
func Nz = N.z;
```

## P
Current point.
```ffpp
real cx = P.x;
real cy = P.y;
real cz = P.z;
```

## pi
Pi = 3.14159.
```ffpp
real Pi = pi;
```
This is a real value.

## sparsesolver

Sparse matrix solver.

## verbosity
Verbosity level.
```ffpp
int Verbosity = verbosity;
verbosity = 0;
```
0 = nothing, 1 = few informations, 10 = a lot of informations, ...

This is an integer value.

## version
FreeFem++ version.
```ffpp
cout << version << endl;
```

## x
The $x$ coordinate at the current point.
```ffpp
real CurrentX = x;
```
This is a real value.

## y
The $y$ coordinate at the current point.
```ffpp
real CurrentY = y;
```
This is a real value.

## z
The $z$ coordinate at the current point.
```ffpp
real CurrentY = y;
```
This is a real value.



