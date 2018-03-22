## The output solution formats .sol and .solb

With the keyword savesol, we can store a scalar functions, a scalar FE functions,
a vector fields, a vector FE fields, a symmetric tensor and a symmetric FE tensor..
Such format is used in medit.

**Extension file `:::freefem .sol`**
The first two lines of the file are :

* `:::freefem MeshVersionFormatted 0`

* `:::freefem Dimension` (I) dim

The following fields begin with one of the following keyword:
`:::freefem SolAtVertices`, `:::freefem SolAtEdges`, `:::freefem SolAtTriangles`, `:::freefem SolAtQuadrilaterals`,
`:::freefem SolAtTetrahedra`, `:::freefem SolAtPentahedra`, `:::freefem SolAtHexahedra`.

In each field, we give then in the next line the number of elements in the solutions
(`:::freefem SolAtVertices`: number of vertices, `:::freefem SolAtTriangles`: number of triangles, ...). In other lines, we give the number of solutions, the type of solution (1: scalar, 2: vector, 3: symmetric tensor). And finally, we give the values of the solutions on the elements.

The file must be ended with the keyword End.

The real element of symmetric tensor :

\begin{eqnarray}
\label{savesol.def.symtensor}
ST^{3d}=\left(
\begin{array}{ccc}
ST_{xx}^{3d} & ST_{xy}^{3d} & ST_{xz}^{3d}\\
ST_{yx}^{3d} & ST_{yy}^{3d} & ST_{yz}^{3d} \\
ST_{zx}^{3d} & ST_{zy}^{3d} & ST_{zz}^{3d}
\end{array}
\right)
\qquad \qquad
ST^{2d}= \left(
\begin{array}{cc}
ST_{xx}^{2d} & ST_{xy}^{2d} \\
ST_{yx}^{2d} & ST_{yy}^{2d}
\end{array}
\right)
\end{eqnarray}

stored in the extension `:::freefem .sol` are respectively $ST_{xx}^{3d}, ST_{yx}^{3d}, ST_{yy}^{3d}, ST_{zx}^{3d}, ST_{zy}^{3d}, ST_{zz}^{3d}$
and  $ST_{xx}^{2d}, ST_{yx}^{2d}, ST_{yy}^{2d}$

An example of field with the keyword `:::freefem SolAtTetrahedra`:

* `:::freefem SolAtTetrahedra`

	(I) NbOfTetrahedrons

	$ \mathtt{ \quad nbsol \quad typesol^1 \quad ... \quad typesol^n }  $
	$\left(\left(\left( \mathtt{U}_{ij}^k, \quad \forall i \in \{1,...,\mathtt{nbrealsol}^k\}\right), %
\quad \forall k \in \{1,...\mathtt{nbsol}\}\right) %
 \quad \forall j \in \{1,...,\mathtt{NbOfTetrahedrons}\}\right)$

where

* $\mathtt{nbsol}$ is an integer equal to the number of solutions

* $\mathtt{typesol^k}$, type of the solution  number $k$, is
	* $\mathtt{typesol^k = 1}$ the solution k is scalar.
	* $\mathtt{typesol^k = 2}$ the solution k is vectorial.
	* $\mathtt{typesol^k = 3}$ the solution k is a symmetric tensor or symmetric matrix.

* $\mathtt{nbrealsol^k}$ number of real to describe solution number $k$ is
	* $\mathtt{nbrealsol^k = 1}$ the solution k is scalar.
	* $\mathtt{nbrealsol^k = dim}$ the solution k is vectorial ($dim$ is the dimension of the solution).
	* $\mathtt{nbrealsol^k = dim*(dim+1)/2}$ the solution k is a symmetric tensor or symmetric matrix.

* $U_{ij}^k$ is a real equal to the value of the component $i$ of the solution $k$ at tetrahedra $j$ on the associated mesh.


This field is written with the notation of Section \ref{meshformatfile.mesh} $\codered$.
The format `:::freefem .solb` is the same as format `:::freefem .sol` but in binary (read/write is faster, storage is less).

A real scalar functions $f1$, a vector fields $\Phi=[\Phi1,\Phi2,\Phi3]$ and a symmetric tensor $ST^{3d}$
(\ref{savesol.def.symtensor}) $\codered$ at the vertices of the three dimensional mesh `:::freefem Th3` is stored in the file `:::freefem f1PhiTh3.sol` using :

```freefem
savesol("f1PhiST3dTh3.sol",Th3, $f1$, [Phi(1), Phi(2), Phi(3), VV3, order=1);
```

where $VV3=[ST_{xx}^{3d}, ST_{yx}^{3d}, ST_{yy}^{3d}, ST_{zx}^{3d}, ST_{zy}^{3d}, ST_{zz}^{3d}]$.
For a two dimensional mesh Th, A real scalar functions $f2$, a vector fields $\Psi=[\Psi1,\Psi2]$ and a symmetric tensor $ST^{2d}$
(\ref{savesol.def.symtensor}) $\codered$ at triangles is stored in the file `:::freefem f2PsiST2dTh3.solb` using :

```freefem
savesol("f2PsiST2dTh3.solb",Th, f2, [Psi(1), Psi(2)], VV2, order=0);
```

where $VV2=[ST_{xx}^{2d}, ST_{yx}^{2d}, ST_{yy}^{2d}]$
The arguments of `:::freefem savesol` functions are the name of a file, a mesh and solutions. These arguments must be given in this order.

The parameters of this keyword are :

* `:::freefem order =` 0 is the solution is given at the center of gravity of elements. 1 is the solution is given at the vertices of elements.

In the file, solutions are stored in this order : scalar solutions, vector solutions and finally symmetric tensor solutions.
