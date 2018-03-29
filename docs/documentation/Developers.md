## File formats

<!-- ### Mesh files -->
<!-- \section{\setS{Mesh Files}}
 \def\Chars#1{{\tt (C*)}  #1}
 \def\Char#1{{\tt (C)}  #1}
 \def\Int#1{ {\tt(I)} #1}
 \def\Real#1{{\tt(R)} #1}
 \def\Bool#1{{\tt(B)} #1}
 \def\Vertex#1{{{\tt @@Vertex}#1}}
 \def\Edge#1{{{\tt @@Edge}#1}}
 \def\Triangle#1{{{\tt @@Tria}#1}}
 \def\Quadrangle#1{{{\tt @@Quad}#1}}
 \def\Tetrahedron#1{{{\tt @@Tetra}#1}}
 \def\Hexahedron#1{{{\tt @@Hexa}#1}}
 \def\Pentahedron#1{{{\tt @@Penta}#1}}
 \def\Loop#1#2{{\bf\Large(}\,#1\,{\bf\Large{,\,\,}}\,#2\,{\bf\Large)}}
 \def\requis{\hfill {\it  requis}}
 \def\facultatif{\quad\quad facultatif}
 \def\need#1{\hfill{\it  requiert le champ\,:\,#1}} -->

### Mesh file data structure
The mesh data structure, output of a mesh generation algorithm, refers to the geometric data structure and in some case to another
mesh data structure.

In this case, the fields are

```cpp
MeshVersionFormatted 0

Dimension [DIM](int)

Vertices
[Number of vertices](int)
X_1(double) Y_1(double) (Z_1(double)) Ref_1(int)
...
X_nv(double) Y_nv(double) (Z_nv(double)) Ref_nv(int)

Edges
[Number of edges](int)
Vertex1_1(int) Vertex2_1(int) Ref_1(int)
...
Vertex1_ne(int) Vertex2_ne(int) Ref_ne(int)

Triangles
[Number of triangles](int)
Vertex1_1(int) Vertex2_1(int) Vertex3_1(int) Ref_1(int)
...
Vertex1_nt(int) Vertex2_nt(int) Vertex3_nt(int) Ref_nt(int)

Quadrilaterals
[Number of Quadrilaterals](int)
Vertex1_1(int) Vertex2_1(int) Vertex3_1(int) Vertex4_1(int) Ref_1(int)
...
Vertex1_nq(int) Vertex2_nq(int) Vertex3_nq(int) Vertex4_nq(int) Ref_nq(int)

Geometry
[File name of geometric support](char*)

	VertexOnGeometricVertex
	[Number of vertex on geometric vertex](int)
	Vertex_1(int) VertexGeometry_1(int)
	...
	Vertex_nvg(int) VertexGeometry_nvg(int)

	EdgeOnGeometricEdge
	[Number of geometric edge](int)
	Edge_1(int) EdgeGeometry_1(int)
	...
	Edge_neg(int) EdgeGeometry_neg(int)

CrackedEdges
[Number of cracked edges](int)
Edge1_1(int) Edge2_1(int)
...
Edge1_nce(int) Edge2_nce(int)
```

When the current mesh refers to a previous mesh, we have in addition

```cpp
MeshSupportOfVertices
[File name of mesh support](char*)

	VertexOnSupportVertex
	[Number of vertex on support vertex](int)
	Vertex_1(int) VertexSupport_1(int)
	...
	Vertex_nvsv(int) VertexSupport_nvsv(int)

	VertexOnSupportEdge
	[Number of vertex on support edge](int)
	Vertex_1(int) EdgeSupport_1(int) USupport_1(double)
	...
	Vertex_nvse(int) EdgeSupport_nvse(int) USupport_nvse(double)

	VertexOnSupportTriangle
	[Number of vertex on support triangle](int)
	Vertex_1(int) TriangleSupport_1(int) USupport_1(double) VSupport_1(double)
	...
	Vertex_nvst(int) TriangleSupport_nvst(int) USupport_nvst(double) VSupport_nvst(double)

	VertexOnSupportQuadrilaterals
	[Number of vertex on support quadrilaterals]
	Vertex_1(int) TriangleSupport_1(int) USupport_1(double) VSupport_1(double)
	...
	Vertex_nvsq(int) TriangleSupport_nvsq(int) USupport_nvsq(double) VSupport_nvsq(double)
```

 * `nv` means the number of vertices
 * `ne` means the number of edges
 * `nt` means the number of triangles
 * `nq` means the number of quadrilaterals
 * `nvg` means the number of vertex on geometric vertex
 * `neg` means the number of edges on geometric edge
 * `nce` means the number of cracked edges

### `bb` file type to Store Solutions

The file is formatted such that:

```cpp
2 [Number of solutions](int) [Number of vertices](int) 2

U_1_1(double) ... U_ns_1(double)
...
U_1_nv(double) ... U_ns_nv(double)
```

 * `ns` means the number of solutions
 * `nv` means the number of vertices
 * `U_i_j` is the solution component `i` at the vertex `j` on the associated mesh.


### `BB` file type to store solutions

The file is formatted such that:

```cpp
2 [Number of solutions](int) [Type 1](int) ... [Type ns](int) [Number of vertices](int) 2

U_1_1_1(double) ... U_(type_k)_1_1(double)
...
U_1_1_1(double) ... U_(type_k)_nbv_1(double)

...

U_1_1_ns(double) ... U_(type_k)_1_ns(double)
...
U_1_nbv_ns(double) ... U_(type_k)_nbv_ns(double)
```

 * `ns` means the number of solutions
 * `type_k` mean the type of solution `k`:
	 - 1: the solution is scalar (1 value per vertex)
	 - 2: the solution is verctorial (2 value sper vertex)
	 - 3: the solution is a $2\times 2$ symmetric matrix (3 vallues per vertex)
	 - 4: the solution is a $2\times 2$ matrix (4 values per vertex)
 * `nbv` means the number of vertices
 * `U_i_j_k` is the value of the component `i`of the solution `k` at vertex `j` on the associated mesh

<!--
\subsection{Metric File}
 A metric file can be of two types, isotropic or anisotropic.
\label{Metric file}

the isotropic file is such that
{\tt \obeylines
   nbv  1
   h$_i \quad \forall i \in \{1,...,\mathtt{nbv}\}$
}


where
\begin{itemize}
\item {\tt  nbv} is  a integer equal to the number of vertices.
\item   {\tt h$_i$} is the wanted mesh size near the vertex $i$ on background mesh,
the metric is $\mathcal{M}_i=h_i^{-2} Id$, where $ Id $ is the identity matrix.
\end{itemize}

The metric anisotrope
{\tt \obeylines
   nbv  3
   a11$_i$,a21$_i$,a22$_i \quad \forall i \in \{1,...,\mathtt{nbv}\}$
}


where
\begin{itemize}
\item   {\tt nbv} is  a integer equal to the number of vertices,
\item  a11$_i$, a12$_i$, a22$_i$ is metric
$\mathcal{M}_i = \left(\begin{smallmatrix} a11_i & a12_i \\ a12_i & a22_i \end{smallmatrix}\right)$ which define the wanted mesh size
in a vicinity of  the vertex $i$
such that $h$ in direction $u \in \R^2$ is equal to $ |u|/\sqrt{u\cdot\mathcal{M}_i\, u}$ , where $\cdot$ is the dot product
in $\R^2$, and $|\cdot|$ is the classical norm.

\end{itemize}

\subsection{List of  AM\_FMT, AMDBA Meshes}
 \index{file!am}\index{file!am\_fmt}\index{file!amdba}
 The mesh is only composed of triangles and can be defined with the help of
the following two integers and four arrays:

  \begin{ttlist}
  \item [nbt] is the number of triangles.
  \item [nbv] is the number of vertices.

  \item [nu(1:3,1:nbt)] is an integer array giving the three vertex numbers

counterclockwise for each triangle.

  \item [c(1:2,nbv)]    is a real array giving the two coordinates of each vertex.
  \item [refs(nbv)]     is an integer array giving the reference numbers of the
vertices.
  \item [reft(nbv)]     is an integer array giving the reference numbers of the
triangles.
  \end{ttlist}

\paragraph{AM\_FMT Files}\label{AMFMT}
\index{file!am\_fmt}
In fortran the  {\tt am\_fmt}  files are read as follows:

\begin{verbatim}
     open(1,file='xxx.am_fmt',form='formatted',status='old')
       read (1,*) nbv,nbt
       read (1,*)  ((nu(i,j),i=1,3),j=1,nbt)
       read (1,*)  ((c(i,j),i=1,2),j=1,nbv)
       read (1,*)  ( reft(i),i=1,nbt)
       read (1,*)  ( refs(i),i=1,nbv)
     close(1)
\end{verbatim}

\paragraph{AM Files}\label{AM}
\index{file!am}
In fortran the  {\tt am}  files are read as follows:

\begin{verbatim}
     open(1,file='xxx.am',form='unformatted',status='old')
       read (1,*) nbv,nbt
       read (1)  ((nu(i,j),i=1,3),j=1,nbt),
     &   ((c(i,j),i=1,2),j=1,nbv),
     &   ( reft(i),i=1,nbt),
     &   ( refs(i),i=1,nbv)
     close(1)
\end{verbatim}
\paragraph{AMDBA Files}\label{AMDBA}
\index{file!amdba}
In fortran the  {\tt amdba}  files are read as follows:
\begin{verbatim}
     open(1,file='xxx.amdba',form='formatted',status='old')
       read (1,*) nbv,nbt
       read (1,*) (k,(c(i,k),i=1,2),refs(k),j=1,nbv)
       read (1,*) (k,(nu(i,k),i=1,3),reft(k),j=1,nbt)
     close(1)
\end{verbatim}
\paragraph{msh Files}\label{MSH}
\index{file!msh}
First, we add the notions of boundary edges
  \begin{ttlist}
  \item [nbbe] is the number of boundary edge.
  \item [nube(1:2,1:nbbe)] is an integer array giving the two vertex numbers
  \item [refbe(1:nbbe)] is an integer array giving the two vertex numbers
  \end{ttlist}
In fortran the  {\tt msh}  files are read as follows:
\begin{verbatim}
     open(1,file='xxx.msh',form='formatted',status='old')
       read (1,*) nbv,nbt,nbbe
       read (1,*) ((c(i,k),i=1,2),refs(k),j=1,nbv)
       read (1,*) ((nu(i,k),i=1,3),reft(k),j=1,nbt)
       read (1,*) ((ne(i,k),i=1,2), refbe(k),j=1,nbbe)
     close(1)
\end{verbatim}
\paragraph{ftq Files}\label{FTQ}
\index{file!ftq}
In fortran the  {\tt ftq}  files are read as follows:
\begin{verbatim}
     open(1,file='xxx.ftq',form='formatted',status='old')
      read (1,*) nbv,nbe,nbt,nbq
      read (1,*) (k(j),(nu(i,j),i=1,k(j)),reft(j),j=1,nbe)
      read (1,*) ((c(i,k),i=1,2),refs(k),j=1,nbv)
     close(1)
\end{verbatim}
where   if {\tt  k(j) = 3} then the element $j$  is  a triangle and if {\tt k = 4}
the the element $j$   is a quadrilateral.

### The output solution formats .sol and .solb

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

In the file, solutions are stored in this order : scalar solutions, vector solutions and finally symmetric tensor solutions. -->
