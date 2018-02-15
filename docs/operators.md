## +
Addition operator.
```ffpp
real a = 1. + 2.;
```
Works for int, real, complex, string, mesh, mesh3, array.

## -
Substraction operator.
```ffpp
real a = 1. - 2.;
```
Works for int, real, complex, array.

## *
Multiplication operator.
```ffpp
real[int] b;
matrix A
real[int] x = A^-1*b;
```
Works for int, real, complex, array, matrix.

## /
Division operator.
```ffpp
real a = 1. / 2.;
```
Works for int, real, complex.

## ^
Power operator.
```ffpp
real a = 2.^2;
```
Works for int, real, complex, matrix.

In the case of matrix, `^-1` mean for the inverse matrix.

## '
Transpose operator.
```ffpp
real[int] a = b';
```
Works for array and matrix.

## :
Tensor scalar product.
$$
A:B = \sum_{i,j}{A_{ij}B_{ij}}
$$



