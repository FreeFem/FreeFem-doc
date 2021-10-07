.. role:: freefem(code)
   :language: freefem

.. role:: cpp(code)
  :language: cpp

Developers
==========

File formats
------------

.. _meshFileDataStructure:

Mesh file data structure
~~~~~~~~~~~~~~~~~~~~~~~~

The mesh data structure, output of a mesh generation algorithm, refers to the geometric data structure and in some case to another mesh data structure.

In this case, the fields are

.. code-block:: cpp
   :linenos:

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

When the current mesh refers to a previous mesh, we have in addition

.. code-block:: cpp
   :linenos:

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

-  ``nv`` means the number of vertices
-  ``ne`` means the number of edges
-  ``nt`` means the number of triangles
-  ``nq`` means the number of quadrilaterals
-  ``nvg`` means the number of vertex on geometric vertex
-  ``neg`` means the number of edges on geometric edge
-  ``nce`` means the number of cracked edges

bb file type to Store Solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The file is formatted such that:

.. code-block:: cpp
   :linenos:

   2 [Number of solutions](int) [Number of vertices](int) 2

   U_1_1(double) ... U_ns_1(double)
   ...
   U_1_nv(double) ... U_ns_nv(double)

-  ``ns`` means the number of solutions
-  ``nv`` means the number of vertices
-  ``U_i_j`` is the solution component ``i`` at the vertex ``j`` on the associated mesh.

BB file type to store solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The file is formatted such that:

.. code-block:: cpp
   :linenos:

   2 [Number of solutions](int) [Type 1](int) ... [Type ns](int) [Number of vertices](int) 2

   U_1_1_1(double) ... U_(type_k)_1_1(double)
   ...
   U_1_1_1(double) ... U_(type_k)_nbv_1(double)

   ...

   U_1_1_ns(double) ... U_(type_k)_1_ns(double)
   ...
   U_1_nbv_ns(double) ... U_(type_k)_nbv_ns(double)

-  ``ns`` means the number of solutions
-  ``type_k`` mean the type of solution ``k``:

   -  1: the solution is scalar (1 value per vertex)
   -  2: the solution is vectorial (2 values per vertex)
   -  3: the solution is a :math:`2\times 2` symmetric matrix (3 values per vertex)
   -  4: the solution is a :math:`2\times 2` matrix (4 values per vertex)

-  ``nbv`` means the number of vertices
-  ``U_i_j_k`` is the value of the component ``i`` of the solution
   ``k`` at vertex ``j`` on the associated mesh

Metric file
~~~~~~~~~~~

A metric file can be of two types, isotropic or anisotropic.

The isotropic file is such that

.. code-block:: cpp
   :linenos:

   [Number of vertices](int) 1
   h_0(double)
   ...
   h_nv(double)

-  ``nv`` is the number of vertices
-  ``h_i`` is the wanted mesh size near the vertex ``i`` on associated mesh.

The metric is :math:`\mathcal{M}_i = h_i^{-2}I` where :math:`I` is the identity matrix.

The anisotropic file is such that

.. code-block:: cpp
   :linenos:

   [Number of vertices](int) 3
   a11_0(double) a21_0(double) a22_0(double)
   ...
   a11_nv(double) a21_nv(double) a22_nv(double)

-  ``nv`` is the number of vertices
-  ``a11_i``, ``a21_i`` and ``a22_i`` represent metric :math:`\mathcal{M}_i = \left(\begin{array}{cc}a_{11,i} & a_{12,i}\\a{12}_i & a_{22,i}\end{array}\right)` which define the wanted size in a vicinity of the vertex ``i`` such that :math:`h` in direction :math:`u \in \mathbb{R}^2` is equal to :math:`|u|/\sqrt{u\cdot\mathcal{M}_i\, u}`, where :math:`\cdot` is the dot product in :math:`\mathbb{R}^2`, and :math:`|\cdot|` is the classical norm.

List of AM_FMT, AMDBA Meshes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The mesh is only composed of triangles and can be defined with the help of the following two integers and four arrays:

-  ``nbt`` the number of triangles
-  ``nbv`` the number of vertices
-  ``nu(1:3, 1:nbt)`` an integer array giving the three vertex numbers counterclockwise for each triangle
-  ``c(1:2, 1:nbv)`` a real array giving tje two coordinates of each vertex
-  ``refs(1:nbv)`` an integer array giving the reference numbers of the vertices
-  ``reft(1:nbt)`` an integer array giving the reference numbers of the triangles

**AM_FMT Files**

In ``Fortran`` the ``am_fmt`` files are read as follows:

.. code-block:: fortran
   :linenos:

   open (1, file='xxx.am_fmt', form='formatted', status='old')
   read (1, *) nbv, nbt
   read (1, *) ((nu(i, j), i=1, 3), j=1, nbt)
   read (1, *) ((c(i, j), i=1, 2), j=1, nbv)
   read (1, *) ( reft(i), i=1, nbt)
   read (1, *) ( refs(i), i=1, nbv)
   close(1)

**AM Files**

In ``Fortran`` the ``am`` files are read as follows:

.. code-block:: fortran
   :linenos:

   open (1, file='xxx.am', form='unformatted', status='old')
   read (1, *) nbv, nbt
   read (1) ((nu(i, j), i=1, 3), j=1, nbt),
   & ((c(i, j), i=1, 2), j=1, nbv),
   & (reft(i), i=1, nbt),
   & (refs(i), i=1, nbv)
   close(1)

**AMDBA Files**

In ``Fortran`` the ``amdba`` files are read as follows:

.. code-block:: fortran
   :linenos:

   open (1, file='xxx.amdba', form='formatted', status='old')
   read (1, *) nbv, nbt
   read (1, *) (k, (c(i, k), i=1, 2), refs(k), j=1, nbv)
   read (1, *) (k, (nu(i, k), i=1, 3), reft(k), j=1, nbt)
   close(1)

**msh Files**

First, we add the notions of boundary edges

-  ``nbbe`` the number of boundary edge
-  ``nube(1:2, 1:nbbe)`` an integer array giving the two vertex numbers of boundary edges
-  ``refbe(1:nbbe)`` an integer array giving the reference numbers of boundary edges

In ``Fortran`` the ``msh`` files are read as follows:

.. code-block:: fortran
   :linenos:

   open (1, file='xxx.msh', form='formatted', status='old')
   read (1, *) nbv, nbt, nbbe
   read (1, *) ((c(i, k), i=1, 2), refs(k), j=1, nbv)
   read (1, *) ((nu(i, k), i=1, 3), reft(k), j=1, nbt)
   read (1, *) ((ne(i, k), i=1, 2), refbe(k), j=1, nbbe)
   close(1)

**ftq Files**

In ``Fortran`` the ``ftq`` files are read as follows:

.. code-block:: fortran
   :linenos:

   open(1,file='xxx.ftq',form='formatted',status='old')
   read (1,*) nbv,nbe,nbt,nbq
   read (1,*) (k(j),(nu(i,j),i=1,k(j)),reft(j),j=1,nbe)
   read (1,*) ((c(i,k),i=1,2),refs(k),j=1,nbv)
   close(1)

where if ``k(j) = 3`` when the element ``j`` is a triangle and ``k(j) = 4`` when the the element ``j`` is a quadrilateral.

sol and solb files
~~~~~~~~~~~~~~~~~~

With the keyword :freefem:`savesol`, we can store a scalar functions, a scalar finite element functions, a vector fields, a vector finite element fields, a symmetric tensor and a symmetric finite element tensor.

Such format is used in :freefem:`medit`.

**Extension file .sol**

The first two lines of the file are :

-  ``MeshVersionFormatted 0``

-  ``Dimension [DIM](int)``

The following fields begin with one of the following keyword:
``SolAtVertices``, ``SolAtEdges``,
``SolAtTriangles``, ``SolAtQuadrilaterals``,
``SolAtTetrahedra``, ``SolAtPentahedra``,
``SolAtHexahedra``.

In each field, we give then in the next line the number of elements in the solutions (:freefem:`SolAtVertices`: number of vertices, :freefem:`SolAtTriangles`: number of triangles, …).
In other lines, we give the number of solutions, the type of solution (1: scalar, 2: vector, 3: symmetric tensor).
And finally, we give the values of the solutions on the elements.

The file must be ended with the keyword End.

The real element of symmetric tensor :

.. math::
   ST^{3d}=\left(
   \begin{array}{ccc}
      ST_{xx}^{3d} & ST_{xy}^{3d} & ST_{xz}^{3d}\\
      ST_{yx}^{3d} & ST_{yy}^{3d} & ST_{yz}^{3d} \\
      ST_{zx}^{3d} & ST_{zy}^{3d} & ST_{zz}^{3d}
   \end{array}
   \right)
   \quad
   ST^{2d}= \left(
   \begin{array}{cc}
      ST_{xx}^{2d} & ST_{xy}^{2d} \\
      ST_{yx}^{2d} & ST_{yy}^{2d}
   \end{array}
   \right)
   :label: savesol.def.symtensor

stored in the extension :freefem:`.sol` are respectively :math:`ST_{xx}^{3d}, ST_{yx}^{3d}, ST_{yy}^{3d}, ST_{zx}^{3d}, ST_{zy}^{3d}, ST_{zz}^{3d}` and :math:`ST_{xx}^{2d}, ST_{yx}^{2d}, ST_{yy}^{2d}`

An example of field with the keyword :freefem:`SolAtTetrahedra`:

.. code-block:: cpp
   :linenos:

   SolAtTetrahedra
   [Number of tetrahedra](int)
   [Number of solutions](int) [Type of solution 1](int) ... [Type of soution nt](int)

   U_1_1_1(double) ... U_nrs_1_1(double)
   ...
   U_1_ns_1(double) ... U_(nrs_k)_ns_1(double)

   ...

   U_1_1_nt(double) ... U_nrs_1_nt(double)
   ...
   U_1_ns_nt(double) ... U_(nrs_k)_ns_nt(double)

-  ``ns`` is the number of solutions
-  ``typesol_k``, type of the solution number ``k``

   -  ``typesol_k = 1`` the solution ``k`` is scalar
   -  ``typesol_k = 2`` the solution ``k`` is vectorial
   -  ``typesol_k = 3`` the solution ``k`` is a symmetric tensor or symmetric matrix

-  ``nrs_k`` is the number of real to describe solution ``k``

   -  ``nrs_k = 1`` if the solution ``k`` is scalar
   -  ``nrs_k = dim`` if the solution ``k`` is vectorial (``dim`` is the dimension of the solution)
   -  ``nrs_k = dim*(dim+1)/2`` if the solution k is a symmetric tensor or symmetric matrix

-  ``U_i_j_^k`` is a real equal to the value of the component ``i`` of the solution ``k`` at tetrahedron ``j`` on the associated mesh

The format :freefem:`.solb` is the same as format :freefem:`.sol` but in binary (read/write is faster, storage is less).

A real scalar functions :math:`f1`, a vector fields :math:`\mathbf{\Phi} = [\Phi1, \Phi2, \Phi3]` and a symmetric tensor :math:`ST^{3d}` :eq:`savesol.def.symtensor` at the vertices of the three dimensional mesh :freefem:`Th3` is stored in the file :freefem:`f1PhiTh3.sol` using :

.. code-block:: freefem
   :linenos:

   savesol("f1PhiST3dTh3.sol", Th3, f1, [Phi(1), Phi(2), Phi(3)], VV3, order=1);

where :math:`VV3 = [ST_{xx}^{3d}, ST_{yx}^{3d}, ST_{yy}^{3d}, ST_{zx}^{3d}, ST_{zy}^{3d}, ST_{zz}^{3d}]`.

For a two dimensional mesh :freefem:`Th`, A real scalar functions :math:`f2`, a vector fields :math:`\mathbf{\Psi} = [\Psi1, \Psi2]` and a symmetric tensor :math:`ST^{2d}` :eq:`savesol.def.symtensor` at triangles is stored in the file :freefem:`f2PsiST2dTh3.solb` using :

.. code-block:: freefem
   :linenos:

   savesol("f2PsiST2dTh3.solb", Th, f2, [Psi(1), Psi(2)], VV2, order=0);

where :math:`VV2 = [ST_{xx}^{2d}, ST_{yx}^{2d}, ST_{yy}^{2d}]`

The arguments of :freefem:`savesol` functions are the name of a file, a mesh and solutions.
These arguments must be given in this order.

The parameters of this keyword are :

-  :freefem:`order =` 0 is the solution is given at the center of gravity of elements.
   1 is the solution is given at the vertices of elements.

In the file, solutions are stored in this order : scalar solutions, vector solutions and finally symmetric tensor solutions.

.. _developersAddingFiniteElement:

Adding a new finite element
---------------------------

Some notations
~~~~~~~~~~~~~~

For a function :math:`\boldsymbol{f}` taking value in :math:`\mathbb{R}^{N},\, N=1,2,\cdots`, we define the finite element approximation :math:`\Pi_h\boldsymbol{f}` of :math:`\boldsymbol{f}`.

Let us denote the number of the degrees of freedom of the finite element by :math:`NbDoF`.
Then the :math:`i`-th base :math:`\boldsymbol{\omega}^{K}_{i}` (:math:`i=0,\cdots,NbDoF-1`) of the finite element space has the :math:`j`-th component :math:`\mathbf{\omega}^{K}_{ij}` for :math:`j=0,\cdots,N-1`.

The operator :math:`\Pi_{h}` is called the interpolator of the finite element.

We have the identity :math:`\boldsymbol{\omega}^{K}_{i} = \Pi_{h} \boldsymbol{\omega}^{K}_{i}`.

Formally, the interpolator :math:`\Pi_{h}` is constructed by the following formula:

.. math::
   \Pi_{h} \boldsymbol{f} = \sum_{k=0}^{\mathtt{kPi}-1} \alpha_k \boldsymbol{f}_{j_{k}}(P_{p_{k}}) \boldsymbol{\omega}^{K}_{i_{k}}
   :label: eq-interpo

where :math:`P_{p}` is a set of :math:`npPi` points,

In the formula :eq:`eq-interpo`, the list :math:`p_{k},\, j_{k},\, i_{k}` depend just on the type of finite element (not on the element), but the coefficient :math:`\alpha_{k}` can be depending on the element.

.. tip:: Classical scalar Lagrange finite element

   With the classical scalar Lagrange finite element, we have :math:`\mathtt{kPi}=\mathtt{npPi}=\mathtt{NbOfNode}` and

   * :math:`P_{p}` is the point of the nodal points.
   * the :math:`\alpha_k=1`, because we take the value of the function at the point :math:`P_{k}`.
   * :math:`p_{k}=k` , :math:`j_{k}=k` because we have one node per function.
   * :math:`j_{k}=0` because :math:`N=1`.

.. tip:: The Raviart-Thomas finite element

   .. math::
      RT0_{h} = \{ \mathbf{v} \in H(div) / \forall K \in
      \mathcal{T}_{h} \quad \mathbf{v}_{|K}(x,y) =
      \vecttwo{\alpha_{K}}{\beta_{K}} + \gamma_{K}\vecttwo{x}{y} \}
      :label: eq:RT0-fe

   The degrees of freedom are the flux through an edge :math:`e` of the mesh, where the flux of the function :math:`\mathbf{f} : \mathbb{R}^2 \longrightarrow \mathbb{R}^2` is :math:`\int_{e} \mathbf{f}.n_{e}`, :math:`n_{e}` is the unit normal of edge :math:`e` (this implies a orientation of all the edges of the mesh, for example we can use the global numbering of the edge vertices and we just go to small to large number).

   To compute this flux, we use a quadrature formula with one point, the middle point of the edge.
   Consider a triangle :math:`T` with three vertices :math:`(\mathbf{a},\mathbf{b},\mathbf{c})`.

   Let denote the vertices numbers by :math:`i_{a},i_{b},i_{c}`, and define the three edge vectors :math:`\mathbf{e}^{0},\mathbf{e}^{1},\mathbf{e}^{2}` by :math:`sgn(i_{b}-i_{c})(\mathbf{b}-\mathbf{c})`, :math:`sgn(i_{c}-i_{a})(\mathbf{c}-\mathbf{a})`, :math:`sgn(i_{a}-i_{b})(\mathbf{a}-\mathbf{b})`.

   The three basis functions are:

   .. math::
       \boldsymbol{\omega}^{K}_{0}= \frac{sgn(i_{b}-i_{c})}{2|T|}(x-a),\quad \boldsymbol{\omega}^{K}_{1}= \frac{sgn(i_{c}-i_{a})}{2|T|}(x-b),\quad \boldsymbol{\omega}^{K}_{2}= \frac{sgn(i_{a}-i_{b})}{2|T|}(x-c),

   where :math:`|T|` is the area of the triangle :math:`T`.

   So we have :math:`N=2`, :math:`\mathtt{kPi}=6; \mathtt{npPi}=3;` and:

    * :math:`P_{p} = \left\{\frac{\mathbf{b}+\mathbf{c}}{2}, \frac{\mathbf{a}+\mathbf{c}}{2}, \frac{\mathbf{b}+\mathbf{a}}{2} \right\}`

    * :math:`\alpha_{0}= - \mathbf{e}^{0}_{2}, \alpha_{1}= \mathbf{e}^{0}_{1}`,
       :math:`\alpha_{2}= - \mathbf{e}^{1}_{2}, \alpha_{3}= \mathbf{e}^{1}_{1}`,
       :math:`\alpha_{4}= - \mathbf{e}^{2}_{2}, \alpha_{5}= \mathbf{e}^{2}_{1}` (effectively, the vector
       :math:`(-\mathbf{e}^{m}_{2}, \mathbf{e}^{m}_{1})` is orthogonal to the edge :math:`\mathbf{e}^{m}= (e^m_{1},e^m_{2})` with
       a length equal to the side of the edge or equal to :math:`\int_{e^m} 1`).
    * :math:`i_{k}=\{0,0,1,1,2,2\}`,
    * :math:`p_{k}=\{0,0,1,1,2,2\}` , :math:`j_{k}=\{0,1,0,1,0,1,0,1\}`.

Which class to add?
~~~~~~~~~~~~~~~~~~~

Add file ``FE_ADD.cpp`` in directory ``FreeFem-sources/src/femlib`` for
example first to initialize :

.. code-block:: cpp
   :linenos:

   #include "error.hpp"
   #include "rgraph.hpp"
   using namespace std;
   #include "RNM.hpp"
   #include "fem.hpp"
   #include "FESpace.hpp"
   #include "AddNewFE.h"

   namespace Fem2D { ... }

Then add a class which derive for ``public TypeOfFE`` like:

.. code-block:: cpp
   :linenos:

   class TypeOfFE_RTortho : public TypeOfFE { public:
       static int Data[]; //some numbers
       TypeOfFE_RTortho():
       TypeOfFE(
           0+3+0,  //nb degree of freedom on element
           2,      //dimension N of vectorial FE (1 if scalar FE)
           Data,   //the array data
           1,      //nb of subdivision for plotting
           1,      //nb of sub finite element (generaly 1)
           6,      //number kPi of coef to build the interpolator
           3,      //number npPi of integration point to build interpolator
           0       //an array to store the coef \alpha_k to build interpolator
           //here this array is no constant so we have
           //to rebuilt for each element
       )
       {
           const R2 Pt[] = {R2(0.5, 0.5), R2(0.0, 0.5), R2(0.5, 0.0) };
           // the set of Point in hat{K}
           for (int p = 0, kk = 0; p < 3; p++){
               P_Pi_h[p] = Pt[p];
               for (int j = 0; j < 2; j++)
                   pij_alpha[kk++] = IPJ(p, p, j);
           }
       } //definition of i_k, p_k, j_k in interpolator

       void FB(const bool *watdd, const Mesh &Th, const Triangle &K,
           const R2 &PHat, RNMK_ &val) const;

       void Pi_h_alpha(const baseFElement &K, KN_<double> &v) const;
   } ;

where the array data is formed with the concatenation of five array of
size ``NbDoF`` and one array of size ``N``.

This array is:

.. code-block:: cpp
   :linenos:

   int TypeOfFE_RTortho::Data[] = {
       //for each df 0, 1, 3:
       3, 4, 5, //the support of the node of the df
       0, 0, 0, //the number of the df on the node
       0, 1, 2, //the node of the df
       0, 0, 0, //the df come from which FE (generally 0)
       0, 1, 2, //which are the df on sub FE
       0, 0
   }; //for each component j=0, N-1 it give the sub FE associated

where the support is a number :math:`0,1,2` for vertex support, :math:`3,4,5` for edge support, and finally :math:`6` for element support.

The function to defined the function
:math:`\boldsymbol{\omega}^{K}_{i}`, this function return the value of all the basics function or this derivatives in array ``val``, computed at point ``Phat`` on the reference triangle corresponding to point ``R2 P=K(Phat);`` on the current triangle ``K``.

The index :math:`i,j,k` of the array :math:`val(i,j,k)` correspond to:

-  :math:`i` is the basic function number on finite element
   :math:`i \in [0,NoF[`
-  :math:`j` is the value of component :math:`j \in [0,N[`
-  :math:`k` is the type of computed value
   :math:`f(P),dx(f)(P), dy(f)(P), ...\ i \in [0,\mathtt{last\_operatortype}[`.

   .. note:: For optimization, this value is computed only if ``whatd[k]`` is true, and the numbering is defined with

      .. code-block:: cpp
         :linenos:

         enum operatortype {
         op_id = 0,
         op_dx = 1, op_dy = 2,
         op_dxx = 3,op_dyy = 4,
         op_dyx = 5,op_dxy = 5,
         op_dz = 6,
         op_dzz = 7,
         op_dzx = 8, op_dxz = 8,
         op_dzy = 9, op_dyz = 9
         };
         const int last_operatortype = 10;

The shape function:

.. code-block:: cpp
   :linenos:

   void TypeOfFE_RTortho::FB(const bool *whatd, const Mesh &Th, const Triangle & K,
       const R2 &PHat,RNMK_ &val) const
   {
       R2 P(K(PHat));
       R2 A(K[0]), B(K[1]), C(K[2]);
       R l0 = 1 - P.x-P.y;
       R l1 = P.x, l2 = P.y;
       assert(val.N() >= 3);
       assert(val.M() == 2);
       val = 0;
       R a = 1./(2*K.area);
       R a0 = K.EdgeOrientation(0) * a;
       R a1 = K.EdgeOrientation(1) * a;
       R a2 = K.EdgeOrientation(2) * a;

       if (whatd[op_id]){ //value of the function
           assert(val.K() > op_id);
           RN_ f0(val('.', 0,0)); //value first component
           RN_ f1(val('.', 1,0)); //value second component
           f1[0] = (P.x - A.x)*a0;
           f0[0] = -(P.y - A.y)*a0;

           f1[1] = (P.x - B.x)*a1;
           f0[1] = -(P.y - B.y)*a1;

           f1[2] = (P.x - C.x)*a2;
           f0[2] = -(P.y - C.y)*a2;
       }

       if (whatd[op_dx]){ //value of the dx of function
           assert(val.K() > op_dx);
           val(0,1,op_dx) = a0;
           val(1,1,op_dx) = a1;
           val(2,1,op_dx) = a2;
       }
       if (whatd[op_dy]){
           assert(val.K() > op_dy);
           val(0,0,op_dy) = -a0;
           val(1,0,op_dy) = -a1;
           val(2,0,op_dy) = -a2;
       }

       for (int i = op_dy; i < last_operatortype; i++)
           if (whatd[op_dx])
               assert(op_dy);
   }

The function to defined the coefficient :math:`\alpha_{k}`:

.. code-block:: cpp
   :linenos:

   void TypeOfFE_RT::Pi_h_alpha(const baseFElement &K, KN_<double> &v) const
   {
       const Triangle &T(K.T);

       for (int i = 0, k = 0; i < 3; i++){
           R2 E(T.Edge(i));
           R signe = T.EdgeOrientation(i) ;
           v[k++] = signe*E.y;
           v[k++] = -signe*E.x;
       }
   }

Now , we just need to add a new key work in **FreeFEM**.

Two way, with static or dynamic link so at the end of the file, we add:

**With dynamic link** it is very simple (see section :ref:`Dynamical link <developersDynamicalLink>`), just add before the end of :cpp:`FEM2d namespace`:

.. code-block:: cpp
   :linenos:

      static TypeOfFE_RTortho The_TypeOfFE_RTortho;
      static AddNewFE("RT0Ortho", The_TypeOfFE_RTortho);
   } //FEM2d namespace

Try with ``./load.link`` command in `examples++-load/ <https://github.com/FreeFem/FreeFem-sources/tree/master/examples%2B%2B-load>`__ and see ``BernardiRaugel.cpp`` or ``Morley.cpp`` new finite element examples.

**Otherwise with static link** (for expert only), add

.. code-block:: cpp
   :linenos:

   //let the 2 globals variables
   static TypeOfFE_RTortho The_TypeOfFE_RTortho;
   //the name in freefem
   static ListOfTFE typefemRTOrtho("RT0Ortho", &The_TypeOfFE_RTortho);

   //link with FreeFEM do not work with static library .a
   //so add a extern name to call in init_static_FE
   //(see end of FESpace.cpp)
   void init_FE_ADD() { };
   //end
   } //FEM2d namespace

To inforce in loading of this new finite element, we have to add the two new lines close to the end of files ``src/femlib/FESpace.cpp`` like:

.. code-block:: cpp
   :linenos:

   //correct problem of static library link with new make file
   void init_static_FE()
   { //list of other FE file.o
       extern void init_FE_P2h() ;
       init_FE_P2h() ;
       extern void init_FE_ADD(); //new line 1
       init_FE_ADD(); //new line 2
   }

and now you have to change the makefile.

First, create a file ``FE_ADD.cpp`` contening all this code, like in file ``src/femlib/Element_P2h.cpp``, after modify the ``Makefile.am`` by adding the name of your file to the variable ``EXTRA_DIST`` like:

.. code-block:: cpp
   :linenos:

   # Makefile using Automake + Autoconf
   # ----------------------------------
   # Id

   # This is not compiled as a separate library because its
   # interconnections with other libraries have not been solved.

   EXTRA_DIST=BamgFreeFem.cpp BamgFreeFem.hpp CGNL.hpp CheckPtr.cpp        \
   ConjuguedGradrientNL.cpp DOperator.hpp Drawing.cpp Element_P2h.cpp      \
   Element_P3.cpp Element_RT.cpp fem3.hpp fem.cpp fem.hpp FESpace.cpp      \
   FESpace.hpp FESpace-v0.cpp FQuadTree.cpp FQuadTree.hpp gibbs.cpp        \
   glutdraw.cpp gmres.hpp MatriceCreuse.hpp MatriceCreuse_tpl.hpp          \
   MeshPoint.hpp mortar.cpp mshptg.cpp QuadratureFormular.cpp              \
   QuadratureFormular.hpp RefCounter.hpp RNM.hpp RNM_opc.hpp RNM_op.hpp    \
   RNM_tpl.hpp     FE_ADD.cpp

and do in the **FreeFEM** root directory

.. code-block:: bash
   :linenos:

   autoreconf
   ./reconfigure
   make

For codewarrior compilation add the file in the project an remove the flag in panal PPC linker FreeFm++ Setting Dead-strip Static Initializition Code Flag.

.. _developersDynamicalLink:

Dynamical link
--------------

Now, it’s possible to add built-in functionnalites in **FreeFEM** under the three environnents Linux, Windows and MacOS X 10.3 or newer.

It is agood idea to first try the example ``load.edp`` in directory `example++-load <https://github.com/FreeFem/FreeFem-sources/tree/master/examples%2B%2B-load>`__.

You will need to install a ``compiler`` (generally ``g++/gcc`` compiler) to compile your function.

-  Windows Install the ``cygwin`` environnent or the ``mingw`` one
-  MacOs Install the developer tools ``Xcode`` on the apple DVD
-  Linux/Unix Install the correct compiler (``gcc`` for instance)

Now, assume that you are in a shell window (a ``cygwin`` window under Windows) in the directory `example++-load <https://github.com/FreeFem/FreeFem-sources/tree/master/examples%2B%2B-load>`__.

.. note:: In the sub directory ``include``, they are all the **FreeFEM** include file to make the link with **FreeFEM**.

.. note:: If you try to load dynamically a file with command :freefem:`load "xxx"`
    - Under Unix (Linux or MacOs), the file ``xxx.so`` will be loaded so it must be either in the search directory of routine ``dlopen`` (see the environment variable `$LD_LIBRARY_PATH`) or in the current directory, and the suffix ``".so"`` or the prefix ``"./"`` is automatically added.

    - Under Windows, the file `xxx.dll` will be loaded so it must be in the `loadLibary` search directory which includes the directory of the application,

**Compilation of your module:**

The script ``ff-c++`` compiles and makes the link with **FreeFEM**, but be careful, the script has no way to known if you try to compile for a pure Windows environment or for a cygwin environment so to build the load module under cygwin you must add the ``-cygwin`` parameter.

A first example myfunction.cpp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following defines a new function call ``myfunction`` with no parameter, but using the :math:`x,y` current value.

.. code-block:: cpp
   :linenos:

   #include <iostream>
   #include <cfloat>
   using namespace std;
   #include "error.hpp"
   #include "AFunction.hpp"
   #include "rgraph.hpp"
   #include "RNM.hpp"
   #include "fem.hpp"
   #include "FESpace.hpp"
   #include "MeshPoint.hpp"

   using namespace Fem2D;
   double myfunction(Stack stack){
       //to get FreeFEM data
       MeshPoint &mp = *MeshPointStack(stack); //the struct to get x, y, normal, value
       double x = mp.P.x; //get the current x value
       double y = mp.P.y; //get the current y value
       //cout << "x = " << x << " y=" << y << endl;
       return sin(x)*cos(y);
   }

Now the Problem is to build the link with **FreeFEM**, to do that we need two classes, one to call the function ``myfunction``.

All **FreeFEM** evaluable expression must be a ``C++`` ``struct``/``class`` which derivate from ``E_F0``.
By default this expression does not depend of the mesh position, but if they derivate from ``E_F0mps`` the expression depends of the mesh position, and for more details see [HECHT2002]_.

.. code-block:: cpp
   :linenos:

   //A class build the link with FreeFEM
   //generaly this class are already in AFunction.hpp
   //but unfortunatly, I have no simple function with no parameter
   //in FreeFEM depending of the mesh
   template<class R>
   class OneOperator0s : public OneOperator {
       //the class to define and evaluate a new function
       //It must devive from E_F0 if it is mesh independent
       //or from E_F0mps if it is mesh dependent
       class E_F0_F :public E_F0mps {
       public:
           typedef R (*func)(Stack stack);
           func f; //the pointeur to the fnction myfunction
           E_F0_F(func ff) : f(ff) {}
           //the operator evaluation in FreeFEM
           AnyType operator()(Stack stack) const {return SetAny<R>(f(stack));}
       };
       typedef R (*func)(Stack);
       func f;
       public:
           //the function which build the FreeFEM byte code
           E_F0 *code(const basicAC_F0 &) const { return new E_F0_F(f); }
           //the constructor to say ff is a function without parameter
           //and returning a R
           OneOperator0s(func ff) : OneOperator(map_type[typeid(R).name()]),f(ff){}
   };

To finish we must add this new function in **FreeFEM** table, to do that include :

.. code-block:: cpp
   :linenos:

    void init(){
        Global.Add("myfunction", "(", new OneOperator0s<double>(myfunction));
    }
    LOADFUNC(init);

It will be called automatically at load module time.

To compile and link, use the ``ff-c++`` script :

.. code-block:: cpp
   :linenos:

   ff-c++ myfunction.cpp
   g++ -c -g -Iinclude myfunction.cpp
   g++ -bundle -undefined dynamic_lookup -g myfunction.o -o ./myfunction.dylib

To try the simple example under Linux or MacOS, do ``FreeFem++-nw load.edp``

The output must be:

.. code-block:: bash
   :linenos:

   -- FreeFem++ v  *.****** (date *** ** *** ****, **:**:** (UTC+0*00))
    Load: lg_fem lg_mesh lg_mesh3 eigenvalue
       1 : // Example of dynamic function load
       2 : // --------------------------------
       3 : // $Id$
       4 :
       5 :  load "myfunction"
       6 : // dumptable(cout);
       7 :  mesh Th=square(5,5);
       8 :  fespace Vh(Th,P1);
       9 :  Vh uh= myfunction(); // warning do not forget ()
      10 :  cout << uh[].min << " " << uh[].max << endl;
      11 :  cout << " test io ( " << endl;
      12 :  testio();
      13 :  cout << " ) end test io .. " << endl; sizestack + 1024 =1416  ( 392 )

     -- Square mesh : nb vertices  =36 ,  nb triangles = 50 ,  nb boundary edges 20
   0 0.841471
    test io (
    test cout 3.14159
    test cout 512
    test cerr 3.14159
    test cerr 512
    ) end test io ..
   times: compile 0.012854s, execution 0.000313s,  mpirank:0
    CodeAlloc : nb ptr  2715,  size :371104 mpirank: 0
   Ok: Normal End

Under Windows, launch **FreeFEM** with the mouse (or ctrl O) on the example.

Example: Discrete Fast Fourier Transform
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This will add FFT to **FreeFEM**, taken from `FFTW <http://www.fftw.org/>`__. To download and install under ``download/include`` just go in ``download/fftw`` and try ``make``.

The 1D dfft (fast discret fourier transform) for a simple array :math:`f` of size :math:`n` is defined by the following formula:

.. math::
   \mathtt{dfft}(f,\varepsilon)_{k} = \sum_{j=0}^{n-1} f_i e^{\varepsilon 2\pi i kj/n}

The 2D DFFT for an array of size :math:`N=n\times m` is:

.. math::
   \mathtt{dfft}(f,m,\varepsilon)_{k+nl} = \sum_{j'=0}^{m-1} \sum_{j=0}^{n-1} f_{i+nj} e^{\varepsilon 2\pi i (kj/n+lj'/m) }

.. note:: The value :math:`n` is given by :math:`size(f)/m`, and the numbering is row-major order.

So the classical discrete DFFT is :math:`\hat{f}=\mathtt{dfft}(f,-1)/\sqrt{n}` and the reverse dFFT :math:`f=\mathtt{dfft}(\hat{f},1)/\sqrt{n}`

.. note:: The 2D Laplace operator is

    .. math::
           f(x,y) = 1/\sqrt{N} \sum_{j'=0}^{m-1} \sum_{j=0}^{n-1} \hat{f}_{i+nj} e^{\varepsilon 2\pi i (x j+ yj') }

    and we have

    .. math::
           f_{k+nl} = f(k/n,l/m)

    So

    .. math::
           \widehat{\Delta f_{kl}} = -( (2\pi)^2 ( (\tilde{k})^2+(\tilde{l})^2)) \widehat{ f_{kl}} \\

    where :math:`\tilde{k} = k` if :math:`k \leq n/2` else :math:`\tilde{k} = k-n` and :math:`\tilde{l} = l` if :math:`l \leq m/2` else :math:`\tilde{l} = l-m`.

   And to have a real function we need all modes to be symmetric around zero, so :math:`n` and :math:`m` must be odd.

**Compile to build a new library**

.. code-block:: bash
   :linenos:

   ff-c++ dfft.cpp ../download/install/lib/libfftw3.a -I../download/install/include
   export MACOSX_DEPLOYMENT_TARGET=10.3
   g++ -c -Iinclude -I../download/install/include dfft.cpp
   g++ -bundle -undefined dynamic_lookup dfft.o -o ./dfft.dylib ../download/install/lib/libfftw3.a

To test, try :ref:`FFT example <exampleFFT>`.

Load Module for Dervieux P0-P1 Finite Volume Method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The associed edp file is
`examples++-load/convect_dervieux.edp <https://github.com/FreeFem/FreeFem-sources/blob/master/examples%2B%2B-load/convect_dervieux.edp>`__.

See
`mat_dervieux.cpp <https://github.com/FreeFem/FreeFem-sources/blob/master/examples%2B%2B-load/mat_dervieux.cpp>`__.

More on Adding a new finite element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First read the :ref:`Adding a new finite element section <developersAddingFiniteElement>`, we add two new finite elements examples in the directory `examples++-load <https://github.com/FreeFem/FreeFem-sources/tree/master/examples%2B%2B-load>`__.

The Bernardi-Raugel Element
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Bernardi-Raugel finite element is meant to solve the Navier Stokes equations in :math:`u,p` formulation; the velocity space :math:`P^{br}_K` is minimal to prove the inf-sup condition with piecewise constant pressure by triangle.

The finite element space :math:`V_h` is

.. math::
   V_h= \{u\in H^1(\Omega)^2 ; \quad \forall K \in T_h, u_{|K} \in P^{br}_K \}

where

.. math::
   P^{br}_K = span \{ \lambda^K_i e_k \}_{i=1,2,3, k= 1,2} \cup \{ \lambda^K_i\lambda^K_{i+1} n^K_{i+2}\}_{i=1,2,3}

with notation :math:`4=1, 5=2` and where :math:`\lambda^K_i` are the barycentric coordinates of the triangle :math:`K`, :math:`(e_k)_{k=1,2}` the canonical basis of :math:`\mathbb{R}^2` and :math:`n^K_k` the outer normal of triangle :math:`K` opposite to vertex :math:`k`.

See
`BernardiRaugel.cpp <https://github.com/FreeFem/FreeFem-sources/blob/master/examples%2B%2B-load/BernardiRaugel.cpp>`__.

A way to check the finite element

.. code-block:: freefem
   :linenos:

   load "BernardiRaugel"

   // Macro
   //a macro the compute numerical derivative
   macro DD(f, hx, hy) ( (f(x1+hx, y1+hy) - f(x1-hx, y1-hy))/(2*(hx+hy)) ) //

   // Mesh
   mesh Th = square(1, 1, [10*(x+y/3), 10*(y-x/3)]);

   // Parameters
   real x1 = 0.7, y1 = 0.9, h = 1e-7;
   int it1 = Th(x1, y1).nuTriangle;

   // Fespace
   fespace Vh(Th, P2BR);
   Vh [a1, a2], [b1, b2], [c1, c2];


   for (int i = 0; i < Vh.ndofK; ++i)
       cout << i << " " << Vh(0,i) << endl;

   for (int i = 0; i < Vh.ndofK; ++i)
   {
       a1[] = 0;
       int j = Vh(it1, i);
       a1[][j] = 1;
       plot([a1, a2], wait=1);
       [b1, b2] = [a1, a2]; //do the interpolation

       c1[] = a1[] - b1[];
       cout << " ---------" << i << " " << c1[].max << " " << c1[].min << endl;
       cout << " a = " << a1[] <<endl;
       cout << " b = " << b1[] <<endl;
       assert(c1[].max < 1e-9 && c1[].min > -1e-9); //check if the interpolation is correct

       // check the derivative and numerical derivative
       cout << " dx(a1)(x1, y1) = " << dx(a1)(x1, y1) << " == " << DD(a1, h, 0) << endl;
       assert( abs(dx(a1)(x1, y1) - DD(a1, h, 0) ) < 1e-5);
       assert( abs(dx(a2)(x1, y1) - DD(a2, h, 0) ) < 1e-5);
       assert( abs(dy(a1)(x1, y1) - DD(a1, 0, h) ) < 1e-5);
       assert( abs(dy(a2)(x1, y1) - DD(a2, 0, h) ) < 1e-5);
   }

A real example using this finite element, just a small modification of the Navier-Stokes P2-P1 example, just the begenning is change to

.. code-block:: freefem
   :linenos:

   load "BernardiRaugel"

   real s0 = clock();
   mesh Th = square(10, 10);
   fespace Vh2(Th, P2BR);
   fespace Vh(Th, P0);
   Vh2 [u1, u2], [up1, up2];
   Vh2 [v1, v2];

And the plot instruction is also changed because the pressure is constant, and we cannot plot isovalues of peacewise constant functions.

The Morley Element
^^^^^^^^^^^^^^^^^^

See the example
`bilapMorley.edp <https://github.com/FreeFem/FreeFem-sources/blob/master/examples%2B%2B-load/bilapMorley.edp>`__.
