## +
Addition operator.
```freefem
real a = 1. + 2.;
```
Works for `:::freefem int`, `:::freefem real`, `:::freefem complex`, `:::freefem string`, `:::freefem mesh`, `:::freefem mesh3`, array.

## -
Substraction operator.
```freefem
real a = 1. - 2.;
```
Works for `:::freefem int`, `:::freefem real`, `:::freefem complex`, array.

## *
Multiplication operator.
```freefem
real[int] b;
matrix A
real[int] x = A^-1*b;
```
Works for `:::freefem int`, `:::freefem real`, `:::freefem complex`, array, `:::freefem matrix`.

## .*
Term by term multiplication.
```freefem
matrix A = B .* C;
```

## /
Division operator.
```freefem
real a = 1. / 2.;
```
Works for `:::freefem int`, `:::freefem real`, `:::freefem complex`.

## ./
Term by term division.
```freefem
matrix A = B ./ C;
```

## %
Remainder from the division.

```freefem
int a = 1 % 2;
```
Works for `:::freefem int`, `:::freefem real`.

## ^
Power operator.
```freefem
real a = 2.^2;
```
Works for `:::freefem int`, `:::freefem real`, `:::freefem complex`, `:::freefem matrix`.

## ^-1
Inverse of a `:::freefem matrix`.

```freefem
real[int] Res = A^-1 * b;
```

!!!warning
	This operator can not be used to directly create a matrix, see [Matrix inversion](/examples/#matrix-inversion).

## '
Transpose operator.
```freefem
real[int] a = b';
```
Works for array and `:::freefem matrix`.

!!!note
	For `:::freefem matrix<complex>`, the `''` operator return the Hermitian tranpose.

## :
Tensor scalar product.
$$
A:B = \sum_{i,j}{A_{ij}B_{ij}}
$$

## ? :
`C++` arithmetical if expression.

`a ? b : c` is equal to `b` if the `a` is true, `c` otherwise.

!!!example "Example with `:::freefem int`"
	```freefem
	int a = 12;
	int b = 5;

	cout << a << " + " << b << " = " << a + b << endl;
	cout << a << " - " << b << " = " << a - b << endl;
	cout << a << " * " << b << " = " << a * b << endl;
	cout << a << " / " << b << " = " << a / b << endl;
	cout << a << " % " << b << " = " << a % b << endl;
	cout << a << " ^ " << b << " = " << a ^ b << endl;
	cout << "( " << a << " < " << b << " ? " << a << " : " << b << ") = " << (a < b ? a : b) << endl;
	```

	The output of this script is:
	```bash
	12 + 5 = 17
	12 - 5 = 7
	12 * 5 = 60
	12 / 5 = 2
	12 % 5 = 2
	12 ^ 5 = 248832
	( 12 < 5 ? 12 : 5) = 5
	```

!!!example "Example with `:::freefem real`"
	```freefem
	real a = qsrt(2.);
	real b = pi;

	cout << a << " + " << b << " = " << a + b << endl;
	cout << a << " - " << b << " = " << a - b << endl;
	cout << a << " * " << b << " = " << a * b << endl;
	cout << a << " / " << b << " = " << a / b << endl;
	cout << a << " % " << b << " = " << a % b << endl;
	cout << a << " ^ " << b << " = " << a ^ b << endl;
	cout << "( " << a << " < " << b << " ? " << a << " : " << b << ") = " << (a < b ? a : b) << endl;
	```

	The output of this script is:
	```bash
	1.41421 + 3.14159 = 4.55581
	1.41421 - 3.14159 = -1.72738
	1.41421 * 3.14159 = 4.44288
	1.41421 / 3.14159 = 0.450158
	1.41421 % 3.14159 = 1
	1.41421 ^ 3.14159 = 2.97069
	( 1.41421 < 3.14159 ? 1.41421 : 3.14159) = 1.41421
	```
